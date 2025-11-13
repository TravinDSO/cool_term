#!/usr/bin/env python3
"""
Terminal dashboard for the space-operations alpha assistant.

This module orchestrates a curses UI with four coordinated panels: an overview,
activity log, control hints, an inner-system map up top, and an outer-system map
anchored in the bottom-right quadrant (which now also tracks transient visitors
like comet 3I/ATLAS). It glues together the
orbital math helpers, timing loop, and rendering utilities so other parts of
the solution can launch the dashboard to monitor simulation state. The
top-right quadrant focuses on inner planetary positioning while the bottom-right
provides the gas-giant perspective, using multi-scale maps so viewers can reason
about the entire solar system at a glance. The lower-left quadrant now streams
headline snippets from multiple live news feeds, rolling them like a terminal
read-out, while the upper-left quadrant presents a sanitized local weather tile that omits any identifying
location metadata, animates bar gauges between API refreshes, and reserves a powered status strip with a scrolling ASCII texture to keep the display feeling energized.
"""

from __future__ import annotations

# Windows compatibility: try standard curses first, fallback to windows-curses
try:
    import curses
except ImportError:
    try:
        import windows_curses as curses
    except ImportError:
        raise ImportError(
            "curses module not found. On Windows, install it with: pip install windows-curses"
        )

import datetime as dt
import math
import sys
from email.utils import parsedate_to_datetime
import random
import textwrap
import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from urllib.parse import urlencode


# Orbital data: mean longitudes referenced to J2000 (degrees) with circular orbits.
PLANETS = [
    ("Mercury", 87.9691, 0.39, 252.250906, "m"),
    ("Venus", 224.701, 0.72, 181.979801, "V"),
    ("Earth", 365.256, 1.00, 100.466457, "E"),
    ("Mars", 686.980, 1.52, 355.433000, "M"),
    ("Jupiter", 4332.589, 5.20, 34.351519, "J"),
    ("Saturn", 10759.22, 9.58, 50.077444, "S"),
    ("Uranus", 30685.4, 19.20, 314.055005, "U"),
    ("Neptune", 60189.0, 30.05, 304.348665, "N"),
]

ADDITIONAL_INNER_OBJECTS: Tuple[Tuple[str, str], ...] = (
    ("Parker Solar Probe", "p"),
    ("Juice Spacecraft", "j"),
)

ADDITIONAL_OUTER_OBJECTS: Tuple[Tuple[str, str], ...] = (
    ("Europa Clipper", "c"),
)

INNER_PLANET_NAMES = {
    "Mercury",
    "Venus",
    "Earth",
    "Mars",
    "3I/ATLAS",
    "Parker Solar Probe",
    "Juice Spacecraft",
}
OUTER_PLANET_NAMES = {
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
    "3I/ATLAS",
    "Europa Clipper",
}

# NASA HORIZONS API planet IDs for fetching real-time orbital data
HORIZONS_PLANET_IDS: Dict[str, str] = {
    "Mercury": "199",
    "Venus": "299",
    "Earth": "399",
    "Mars": "499",
    "Jupiter": "599",
    "Saturn": "699",
    "Uranus": "799",
    "Neptune": "899",
    "Parker Solar Probe": "PARKER SOLAR PROBE",
    "Juice Spacecraft": "JUICE",
    "Europa Clipper": "EUROPA CLIPPER",
}
HORIZONS_API_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
HORIZONS_TIMEOUT = 10
# Enable debug logging for HORIZONS API (set via environment variable)
HORIZONS_DEBUG = os.getenv("HORIZONS_DEBUG", "false").lower() == "true"

GAUSSIAN_GRAVITATIONAL_CONSTANT = 0.01720209895

COMET_3I_ATLAS = (
    "3I/ATLAS",
    -0.2639,  # semi-major axis (AU, hyperbolic => negative)
    6.1396,  # eccentricity
    175.1131,  # inclination (deg)
    322.157,  # longitude of ascending node (deg)
    128.010,  # argument of perihelion (deg)
    2460977.98,  # time of perihelion passage (Julian Date)
    1.3564,  # perihelion distance (AU)
    "*",  # display symbol
)

MAX_ORBIT_RADIUS = max(p[2] for p in PLANETS)
MIN_ORBIT_RADIUS = min(p[2] for p in PLANETS)
MIN_RADIUS_FRACTION = 0.05  # ensures even Mercury sits away from the sun marker
J2000 = dt.datetime(2000, 1, 1, 12, tzinfo=dt.timezone.utc)

# Terminal cells are taller than they are wide; adjust vertical plotting to keep
# orbital paths visually circular. This ratio can be tuned per font/terminal.
CHAR_ASPECT_RATIO = 0.5

DEFAULT_PLANET_COLORS = {
    "Mercury": curses.COLOR_WHITE,
    "Venus": curses.COLOR_YELLOW,
    "Earth": curses.COLOR_GREEN,
    "Mars": curses.COLOR_RED,
    "Jupiter": curses.COLOR_MAGENTA,
    "Saturn": curses.COLOR_CYAN,
    "Uranus": curses.COLOR_BLUE,
    "Neptune": curses.COLOR_WHITE,
    "3I/ATLAS": curses.COLOR_YELLOW,
    "Parker Solar Probe": curses.COLOR_YELLOW,
    "Juice Spacecraft": curses.COLOR_GREEN,
    "Europa Clipper": curses.COLOR_CYAN,
}

COLOR_SUPPORT = False
PLANET_COLOR_PAIRS: Dict[str, int] = {}

# Environment bootstrap: prioritize the most specific dotenv files first so they
# win when multiple directories contain overlapping key definitions. We never
# overwrite values that are already present in the process environment.
ENV_FILENAMES: Tuple[str, ...] = (".env.local", ".env.development", ".env")

# Weather integration is intentionally keyed off environment variables so users
# can opt-in without hard-coding location metadata in the repo.
WEATHER_LAT_ENV = "ALPHA_WEATHER_LAT"
WEATHER_LON_ENV = "ALPHA_WEATHER_LON"
WEATHER_API_URL = "https://api.open-meteo.com/v1/forecast"
# Request the richest hourly and daily payload Open-Meteo offers so the tile can
# stay informative without storing any location metadata beyond lat/lon.
WEATHER_HOURLY_FIELDS: Tuple[str, ...] = (
    "relativehumidity_2m",
    "apparent_temperature",
    "dewpoint_2m",
    "precipitation",
    "rain",
    "snowfall",
    "cloudcover",
    "visibility",
    "surface_pressure",
    "pressure_msl",
)
WEATHER_DAILY_FIELDS: Tuple[str, ...] = (
    "sunrise",
    "sunset",
    "uv_index_max",
    "precipitation_probability_max",
)
WEATHER_PARAMS: Dict[str, str] = {
    "current_weather": "true",
    "hourly": ",".join(WEATHER_HOURLY_FIELDS),
    "daily": ",".join(WEATHER_DAILY_FIELDS),
    "temperature_unit": "celsius",
    "windspeed_unit": "kmh",
    "timezone": "UTC",
    "forecast_days": "1",
}
# Unified update cycle: all components refresh every 5 minutes
SYNC_UPDATE_INTERVAL_SECONDS = 300
WEATHER_REFRESH_SECONDS = SYNC_UPDATE_INTERVAL_SECONDS
WEATHER_RETRY_SECONDS = SYNC_UPDATE_INTERVAL_SECONDS
WEATHER_SPINNER_CHARS: Tuple[str, ...] = ("-", "\\", "|", "/")
WEATHER_BAR_FILL_CHAR = "#"
WEATHER_BAR_EMPTY_CHAR = "."
WEATHER_BAR_MAX_WIDTH = 28
WEATHER_ANIMATION_MODULO = len(WEATHER_SPINNER_CHARS) * 64
POWER_BAR_TEXTURES: Tuple[str, ...] = (".", ":", "-", "=", "+", "*", "#", "%", "@")
POWER_BAR_MIN_HEIGHT = 3
POWER_BAR_MAX_HEIGHT = 6

# Open-Meteo weather code descriptions sourced from their public documentation.
WEATHER_CODE_MAP: Dict[int, str] = {
    0: "Clear sky",
    1: "Mainly clear",
    2: "Partly cloudy",
    3: "Overcast",
    45: "Fog",
    48: "Depositing rime fog",
    51: "Light drizzle",
    53: "Moderate drizzle",
    55: "Dense drizzle",
    56: "Light freezing drizzle",
    57: "Dense freezing drizzle",
    61: "Slight rain",
    63: "Moderate rain",
    65: "Heavy rain",
    66: "Light freezing rain",
    67: "Heavy freezing rain",
    71: "Slight snow",
    73: "Moderate snow",
    75: "Heavy snow",
    77: "Snow grains",
    80: "Light rain showers",
    81: "Moderate rain showers",
    82: "Violent rain showers",
    85: "Light snow showers",
    86: "Heavy snow showers",
    95: "Thunderstorm",
    96: "Thunderstorm with slight hail",
    99: "Thunderstorm with heavy hail",
}

# News ticker timing tuned for readability while keeping the feed feeling alive.
NEWS_REFRESH_SECONDS = SYNC_UPDATE_INTERVAL_SECONDS
NEWS_CHAR_INTERVAL = 0.05
NEWS_HISTORY_LIMIT = 120
NEWS_MAX_ITEMS = 40

# Curated set of high-signal feeds to keep the lower-left panel populated.
NEWS_FEEDS: Tuple[Tuple[str, str], ...] = (
    ("NASA", "https://www.nasa.gov/rss/dyn/breaking_news.rss"),
    ("ESA", "https://www.esa.int/rssfeed/Our_Activities/Space_News"),
    ("SpaceNews", "https://spacenews.com/feed/"),
)


def load_env_files() -> None:
    """Populate os.environ using dotenv-style files discovered up the tree."""
    seen: set[Path] = set()
    start_dir = Path(__file__).resolve().parent
    search_paths = [start_dir, *start_dir.parents]

    for directory in search_paths:
        for filename in ENV_FILENAMES:
            env_path = directory / filename
            if env_path in seen or not env_path.is_file():
                continue
            seen.add(env_path)
            try:
                for raw_line in env_path.read_text(encoding="utf-8").splitlines():
                    stripped = raw_line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    if stripped.lower().startswith("export "):
                        stripped = stripped[7:].lstrip()
                    if "=" not in stripped:
                        continue
                    key, value = stripped.split("=", 1)
                    key = key.strip()
                    value = value.strip()
                    if "#" in value:
                        value = value.split("#", 1)[0].rstrip()
                    if (
                        value.startswith(("'", '"'))
                        and value.endswith(("'", '"'))
                        and len(value) >= 2
                    ):
                        value = value[1:-1]
                    value = value.strip()
                    if key and key not in os.environ:
                        os.environ[key] = value
            except OSError:
                continue


load_env_files()


def wrap_text_lines(text: str, width: int) -> List[str]:
    if width <= 0:
        return [text]
    wrapped = textwrap.wrap(text, width=width, break_long_words=True, break_on_hyphens=False)
    return wrapped or [""]


def load_weather_coordinates() -> Optional[Tuple[float, float]]:
    lat_str = os.getenv(WEATHER_LAT_ENV)
    lon_str = os.getenv(WEATHER_LON_ENV)
    if not lat_str or not lon_str:
        return None
    try:
        return (float(lat_str), float(lon_str))
    except ValueError:
        return None


def wind_direction_from_degrees(degrees: Optional[float]) -> str:
    if degrees is None:
        return "--"
    bearings = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
    index = int(((degrees % 360) + 22.5) // 45) % len(bearings)
    return bearings[index]


def clamp(value: float, lower: float, upper: float) -> float:
    """Clamp a numeric value to the inclusive [lower, upper] interval."""

    return max(lower, min(upper, value))


def build_weather_bar(
    value: float,
    min_value: float,
    max_value: float,
    width: int,
    animation_phase: int,
) -> str:
    """Create an animated ASCII bar representing a weather metric."""

    if max_value <= min_value:
        return WEATHER_BAR_EMPTY_CHAR * max(1, width)

    bar_width = max(1, width)
    normalized = clamp((value - min_value) / (max_value - min_value), 0.0, 1.0)
    filled_cells = int(round(normalized * bar_width))

    bar_chars = [WEATHER_BAR_EMPTY_CHAR] * bar_width
    for i in range(filled_cells):
        bar_chars[i] = WEATHER_BAR_FILL_CHAR

    pointer_index = min(bar_width - 1, max(0, filled_cells - 1))
    spinner_char = WEATHER_SPINNER_CHARS[animation_phase % len(WEATHER_SPINNER_CHARS)]
    bar_chars[pointer_index] = spinner_char

    return "".join(bar_chars)


def make_weather_gauge_line(
    label: str,
    value: Optional[float],
    min_value: float,
    max_value: float,
    display_value: str,
    content_width: int,
    animation_phase: int,
) -> Optional[str]:
    """Compose a simple label/value line without a progress bar."""

    if not isinstance(value, (int, float)):
        return None

    label_prefix = f"{label}:"
    line = f"{label_prefix} {display_value}".rstrip()

    if len(line) > content_width:
        return line[:content_width]

    return line


def update_power_bar_buffer(buffer: str, width: int) -> str:
    """Scroll the power-status bar using random ASCII textures."""

    if width <= 0:
        return ""

    # Pad or trim the buffer to the requested width.
    if len(buffer) < width:
        buffer = (buffer + " " * width)[-width:]
    elif len(buffer) > width:
        buffer = buffer[-width:]

    chars = list(buffer)
    # Shift existing characters left to create motion.
    chars = chars[1:] + [random.choice(POWER_BAR_TEXTURES)]

    # Randomly flicker a handful of cells for texture.
    for _ in range(max(1, width // 12)):
        idx = random.randrange(0, width)
        chars[idx] = random.choice(POWER_BAR_TEXTURES)

    return "".join(chars)


def fetch_weather_data(lat: float, lon: float) -> Optional[Dict[str, Any]]:
    params = WEATHER_PARAMS.copy()
    params["latitude"] = f"{lat:.4f}"
    params["longitude"] = f"{lon:.4f}"
    url = f"{WEATHER_API_URL}?{urlencode(params)}"

    try:
        with urllib.request.urlopen(url, timeout=10) as response:
            payload = response.read().decode("utf-8")
    except (urllib.error.URLError, TimeoutError, ValueError):
        return None

    try:
        data = json.loads(payload)
    except json.JSONDecodeError:
        return None

    current = data.get("current_weather") or {}
    if not current:
        return None

    current_time = current.get("time")
    timestamp: Optional[dt.datetime] = None
    if isinstance(current_time, str):
        try:
            timestamp = dt.datetime.fromisoformat(current_time)
            if timestamp.tzinfo is None:
                timestamp = timestamp.replace(tzinfo=dt.timezone.utc)
        except ValueError:
            timestamp = None

    hourly = data.get("hourly") or {}
    hourly_index: Optional[int] = None
    humidity = None
    apparent = None
    dewpoint_c = None
    precip_mm = None
    rain_mm = None
    snow_mm = None
    cloud_cover = None
    visibility_m = None
    surface_pressure_hpa = None
    pressure_msl_hpa = None

    hourly_times = hourly.get("time")
    if isinstance(hourly_times, list) and current_time:
        try:
            hourly_index = hourly_times.index(current_time)
        except ValueError:
            hourly_index = None

    def _extract_hourly_numeric(key: str) -> Optional[float]:
        if hourly_index is None:
            return None
        series = hourly.get(key)
        if not isinstance(series, list) or hourly_index >= len(series):
            return None
        value = series[hourly_index]
        if isinstance(value, (int, float)):
            return float(value)
        return None

    humidity = _extract_hourly_numeric("relativehumidity_2m")
    apparent = _extract_hourly_numeric("apparent_temperature")
    dewpoint_c = _extract_hourly_numeric("dewpoint_2m")
    precip_mm = _extract_hourly_numeric("precipitation")
    rain_mm = _extract_hourly_numeric("rain")
    snow_mm = _extract_hourly_numeric("snowfall")
    cloud_cover = _extract_hourly_numeric("cloudcover")
    visibility_m = _extract_hourly_numeric("visibility")
    surface_pressure_hpa = _extract_hourly_numeric("surface_pressure")
    pressure_msl_hpa = _extract_hourly_numeric("pressure_msl")

    temperature_c = current.get("temperature")
    windspeed_kmh = current.get("windspeed")
    wind_dir_deg = current.get("winddirection")
    weather_code = current.get("weathercode")
    is_day = current.get("is_day")

    visibility_km = visibility_m / 1000.0 if isinstance(visibility_m, (int, float)) else None
    visibility_miles = visibility_km * 0.621371 if isinstance(visibility_km, (int, float)) else None
    dewpoint_f = dewpoint_c * 9 / 5 + 32 if isinstance(dewpoint_c, (int, float)) else None
    precip_in = precip_mm * 0.0393701 if isinstance(precip_mm, (int, float)) else None
    rain_in = rain_mm * 0.0393701 if isinstance(rain_mm, (int, float)) else None
    snow_in = snow_mm * 0.0393701 if isinstance(snow_mm, (int, float)) else None
    surface_pressure_inhg = (
        surface_pressure_hpa * 0.0295299830714 if isinstance(surface_pressure_hpa, (int, float)) else None
    )
    pressure_msl_inhg = (
        pressure_msl_hpa * 0.0295299830714 if isinstance(pressure_msl_hpa, (int, float)) else None
    )

    daily = data.get("daily") or {}
    daily_times = daily.get("time")
    daily_index: Optional[int] = None
    if isinstance(daily_times, list) and daily_times:
        if timestamp:
            target_date = timestamp.date().isoformat()
            try:
                daily_index = daily_times.index(target_date)
            except ValueError:
                daily_index = 0
        else:
            daily_index = 0

    def _extract_daily_numeric(key: str) -> Optional[float]:
        if daily_index is None:
            return None
        series = daily.get(key)
        if not isinstance(series, list) or daily_index >= len(series):
            return None
        value = series[daily_index]
        if isinstance(value, (int, float)):
            return float(value)
        return None

    def _extract_daily_datetime(key: str) -> Optional[dt.datetime]:
        if daily_index is None:
            return None
        series = daily.get(key)
        if not isinstance(series, list) or daily_index >= len(series):
            return None
        value = series[daily_index]
        if not isinstance(value, str):
            return None
        try:
            parsed = dt.datetime.fromisoformat(value)
        except ValueError:
            return None
        if parsed.tzinfo is None:
            parsed = parsed.replace(tzinfo=dt.timezone.utc)
        else:
            parsed = parsed.astimezone(dt.timezone.utc)
        return parsed

    sunrise_dt = _extract_daily_datetime("sunrise")
    sunset_dt = _extract_daily_datetime("sunset")
    uv_index_max = _extract_daily_numeric("uv_index_max")
    precip_probability_max = _extract_daily_numeric("precipitation_probability_max")

    temperature_f = temperature_c * 9 / 5 + 32 if isinstance(temperature_c, (int, float)) else None
    apparent_c = apparent if isinstance(apparent, (int, float)) else temperature_c
    apparent_f = apparent_c * 9 / 5 + 32 if isinstance(apparent_c, (int, float)) else None
    windspeed_mph = windspeed_kmh * 0.621371 if isinstance(windspeed_kmh, (int, float)) else None

    return {
        "timestamp": timestamp,
        "temperature_c": temperature_c,
        "temperature_f": temperature_f,
        "apparent_c": apparent_c,
        "apparent_f": apparent_f,
        "windspeed_kmh": windspeed_kmh,
        "windspeed_mph": windspeed_mph,
        "wind_direction_deg": wind_dir_deg,
        "wind_direction_cardinal": wind_direction_from_degrees(wind_dir_deg),
        "humidity": humidity,
        "dewpoint_c": dewpoint_c,
        "dewpoint_f": dewpoint_f,
        "cloud_cover": cloud_cover,
        "visibility_km": visibility_km,
        "visibility_miles": visibility_miles,
        "surface_pressure_hpa": surface_pressure_hpa,
        "surface_pressure_inhg": surface_pressure_inhg,
        "pressure_msl_hpa": pressure_msl_hpa,
        "pressure_msl_inhg": pressure_msl_inhg,
        "precipitation_mm": precip_mm,
        "precipitation_in": precip_in,
        "rain_mm": rain_mm,
        "rain_in": rain_in,
        "snow_mm": snow_mm,
        "snow_in": snow_in,
        "uv_index_max": uv_index_max,
        "precipitation_probability_max": precip_probability_max,
        "sunrise": sunrise_dt,
        "sunset": sunset_dt,
        "is_day": bool(is_day) if is_day is not None else None,
        "condition": WEATHER_CODE_MAP.get(weather_code, "Unknown"),
        "weather_code": weather_code,
    }


def init_color_support() -> None:
    """Configure color pairs for supported bodies when terminal colors are available."""

    global COLOR_SUPPORT, PLANET_COLOR_PAIRS

    if not curses.has_colors():
        COLOR_SUPPORT = False
        PLANET_COLOR_PAIRS.clear()
        return

    try:
        curses.start_color()
        curses.use_default_colors()
    except curses.error:
        COLOR_SUPPORT = False
        PLANET_COLOR_PAIRS.clear()
        return

    COLOR_SUPPORT = True
    PLANET_COLOR_PAIRS.clear()

    for index, (body, color) in enumerate(DEFAULT_PLANET_COLORS.items(), start=1):
        try:
            curses.init_pair(index, color, -1)
        except curses.error:
            COLOR_SUPPORT = False
            PLANET_COLOR_PAIRS.clear()
            return
        PLANET_COLOR_PAIRS[body] = index


def _days_since_epoch(moment: dt.datetime) -> float:
    delta = moment - J2000
    return delta.total_seconds() / 86400.0


def fetch_latest_headlines(max_items: int = NEWS_MAX_ITEMS) -> List[str]:
    """Fetch and format recent headlines aggregated from configured RSS feeds."""

    aggregated: List[Tuple[Optional[dt.datetime], str]] = []

    for source_name, feed_url in NEWS_FEEDS:
        try:
            with urllib.request.urlopen(feed_url, timeout=10) as response:
                data = response.read()
        except (urllib.error.URLError, TimeoutError, ValueError):
            continue

        try:
            root = ET.fromstring(data)
        except ET.ParseError:
            continue

        for item in root.findall(".//item"):
            title = (item.findtext("title") or "").strip()
            if not title:
                continue

            pub_date = item.findtext("pubDate")
            published: Optional[dt.datetime]
            timestamp = "----"
            if pub_date:
                try:
                    parsed = parsedate_to_datetime(pub_date)
                    if parsed.tzinfo is None:
                        parsed = parsed.replace(tzinfo=dt.timezone.utc)
                    published = parsed.astimezone(dt.timezone.utc)
                    timestamp = published.strftime("%H:%MZ")
                except (TypeError, ValueError, OverflowError):
                    published = None
                else:
                    aggregated.append((published, f"[{source_name} {timestamp}] {title}"))
                    continue

            aggregated.append((None, f"[{source_name} {timestamp}] {title}"))

    if not aggregated:
        return []

    aggregated.sort(key=lambda entry: entry[0] or dt.datetime.min.replace(tzinfo=dt.timezone.utc))
    if len(aggregated) > max_items:
        aggregated = aggregated[-max_items:]

    random.shuffle(aggregated)

    return [entry[1] for entry in aggregated[:max_items]]


def _datetime_to_julian_date(moment: dt.datetime) -> float:
    """Convert a datetime into Julian Date, matching the J2000 epoch reference."""

    return 2451545.0 + _days_since_epoch(moment)


def _solve_hyperbolic_anomaly(mean_anomaly: float, eccentricity: float) -> float:
    """Newton-Raphson solver for hyperbolic Kepler equation."""

    if math.isclose(mean_anomaly, 0.0):
        return 0.0

    # Use inverse hyperbolic sine as an efficient initial guess
    if eccentricity <= 1.0:
        raise ValueError("Hyperbolic anomaly solver requires e > 1")

    guess = math.asinh(mean_anomaly / eccentricity)
    F = guess
    for _ in range(60):
        sinh_F = math.sinh(F)
        cosh_F = math.cosh(F)
        f_val = eccentricity * sinh_F - F - mean_anomaly
        f_prime = eccentricity * cosh_F - 1.0
        if math.isclose(f_prime, 0.0, abs_tol=1e-12):
            break
        delta = f_val / f_prime
        F -= delta
        if abs(delta) < 1e-10:
            break
    return F


def compute_comet_position(epoch: dt.datetime) -> Tuple[str, str, float, float] | None:
    """Compute heliocentric x/y coordinates for comet 3I/ATLAS at the given epoch."""

    (
        name,
        semi_major_axis,
        eccentricity,
        inclination_deg,
        long_asc_node_deg,
        arg_perihelion_deg,
        perihelion_jd,
        _perihelion_distance,
        symbol,
    ) = COMET_3I_ATLAS

    jd = _datetime_to_julian_date(epoch)
    delta_t = jd - perihelion_jd

    if semi_major_axis >= 0 or eccentricity <= 1:
        return None

    mean_motion = GAUSSIAN_GRAVITATIONAL_CONSTANT / math.sqrt(abs(semi_major_axis) ** 3)
    mean_anomaly = mean_motion * delta_t

    try:
        hyperbolic_anomaly = _solve_hyperbolic_anomaly(mean_anomaly, eccentricity)
    except ValueError:
        return None

    cosh_F = math.cosh(hyperbolic_anomaly)
    sinh_F = math.sinh(hyperbolic_anomaly)

    radius = semi_major_axis * (eccentricity * cosh_F - 1.0)
    radius = abs(radius)

    true_anomaly = 2.0 * math.atan(
        math.sqrt((eccentricity + 1.0) / (eccentricity - 1.0)) * math.tanh(hyperbolic_anomaly / 2.0)
    )

    inclination = math.radians(inclination_deg)
    long_asc_node = math.radians(long_asc_node_deg)
    arg_perihelion = math.radians(arg_perihelion_deg)

    cos_O = math.cos(long_asc_node)
    sin_O = math.sin(long_asc_node)
    cos_i = math.cos(inclination)
    sin_i = math.sin(inclination)

    argument_sum = arg_perihelion + true_anomaly
    cos_arg = math.cos(argument_sum)
    sin_arg = math.sin(argument_sum)

    x = radius * (cos_O * cos_arg - sin_O * sin_arg * cos_i)
    y = radius * (sin_O * cos_arg + cos_O * sin_arg * cos_i)

    return (name, symbol, x, y)


def fetch_horizons_position(planet_name: str, epoch: dt.datetime) -> Optional[Tuple[float, float, float, float]]:
    """
    Fetch real-time heliocentric position and velocity from NASA HORIZONS API.
    Returns (x, y, vx, vy) in AU and AU/day, or None if fetch fails.
    
    Note: HORIZONS API response format may vary. This function attempts to parse
    the standard format, but if the API changes or returns unexpected data,
    the function will return None and trigger fallback to calculated positions.
    """
    if planet_name not in HORIZONS_PLANET_IDS:
        return None
    
    planet_id = HORIZONS_PLANET_IDS[planet_name]
    command_value = str(planet_id).strip()
    if (
        isinstance(command_value, str)
        and not command_value.startswith("\"")
        and " " in command_value
        and not re.fullmatch(r"-?\d+", command_value)
    ):
        command_value = f'"{command_value}"'
    
    # HORIZONS API date format issue: spaces in dates cause parsing errors
    # Try using JD (Julian Date) format or ISO format without spaces
    # Based on error, let's use JD format which is more reliable
    jd = _datetime_to_julian_date(epoch)
    jd_start = f"JD{jd:.6f}"
    # STOP_TIME must be later than START_TIME - add 1 minute (1/1440 of a day)
    jd_stop = f"JD{jd + (1/1440):.6f}"
    
    # HORIZONS API parameters - use JD format for dates to avoid format issues
    # Remove OBJ_DATA to get ephemeris instead of object properties
    # Simplify parameters to ensure ephemeris is generated
    params = {
        "format": "text",
        "COMMAND": command_value,
        "OBJ_DATA": "NO",  # Don't request object data, just ephemeris
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "VECTORS",
        "CENTER": "500@10",  # Sun center (500) at origin (10)
        "START_TIME": jd_start,  # Use Julian Date format
        "STOP_TIME": jd_stop,  # Must be later than START_TIME
        "STEP_SIZE": "1m",
        "VEC_TABLE": "2",  # State vectors (position and velocity)
        "VEC_LABELS": "YES",
        "CSV_FORMAT": "NO",
        # Remove QUANTITIES - let it use defaults for vectors
    }
    
    try:
        url = f"{HORIZONS_API_URL}?{urlencode(params)}"
        
        if HORIZONS_DEBUG:
            # Debug: print the URL being requested (first planet only to avoid spam)
            if planet_name == "Mercury":
                print(f"DEBUG: HORIZONS URL: {url[:200]}...", file=sys.stderr)
        
        with urllib.request.urlopen(url, timeout=HORIZONS_TIMEOUT) as response:
            result_text = response.read().decode("utf-8")
        
        if HORIZONS_DEBUG and planet_name == "Mercury":
            # Debug: save first 500 chars of response
            print(f"DEBUG: HORIZONS Response (first 500 chars):\n{result_text[:500]}", file=sys.stderr)
        
        # Check for API errors in response
        if "API Error" in result_text or "No matches found" in result_text or "Unknown" in result_text:
            if HORIZONS_DEBUG and planet_name == "Mercury":
                print(f"DEBUG: API Error detected in response", file=sys.stderr)
            return None
        
        # Parse HORIZONS response
        # Find the data section (starts after "$$SOE" and ends at "$$EOE")
        # Format example:
        # 2460992.390619000 = A.D. 2025-Nov-12 21:22:29.4816 TDB
        #  X = 9.439702498868449E+07 Y = 1.140655766260261E+08 Z =-7.432645419202745E+03
        #  VX=-2.342648150829879E+01 VY= 1.888918290715260E+01 VZ= 1.197204814396002E-05
        # Note: X, Y, Z are on one line, VX, VY, VZ are on the NEXT line
        # Coordinates are in KM, velocities in km/s, need to convert to AU and AU/day
        KM_TO_AU = 1.0 / 149597870.7
        
        lines = result_text.split("\n")
        in_data = False
        for i, line in enumerate(lines):
            if "$$SOE" in line:
                in_data = True
                continue
            if "$$EOE" in line:
                break
            if in_data and line.strip():
                # Look for the line with X =, Y =, Z = pattern
                # This is the actual data line (not the JD timestamp line)
                line_upper = line.upper()
                if ("X =" in line_upper or "X=" in line_upper) and i + 1 < len(lines):
                    # Extract X and Y from current line
                    x_match = re.search(r'X\s*=\s*([+-]?\d+\.?\d*[Ee][+-]?\d+|[+-]?\d+\.?\d*)', line_upper)
                    y_match = re.search(r'Y\s*=\s*([+-]?\d+\.?\d*[Ee][+-]?\d+|[+-]?\d+\.?\d*)', line_upper)
                    
                    if x_match and y_match:
                        # Get the next line for velocities
                        next_line_upper = lines[i + 1].upper() if i + 1 < len(lines) else ""
                        vx_match = re.search(r'VX\s*=\s*([+-]?\d+\.?\d*[Ee][+-]?\d+|[+-]?\d+\.?\d*)', next_line_upper)
                        vy_match = re.search(r'VY\s*=\s*([+-]?\d+\.?\d*[Ee][+-]?\d+|[+-]?\d+\.?\d*)', next_line_upper)
                        
                        if vx_match and vy_match:
                            try:
                                # Values are in kilometers and km/s, convert to AU and AU/day
                                x_km = float(x_match.group(1))
                                y_km = float(y_match.group(1))
                                vx_kms = float(vx_match.group(1))
                                vy_kms = float(vy_match.group(1))
                                
                                # Convert to AU and AU/day
                                x_au = x_km * KM_TO_AU
                                y_au = y_km * KM_TO_AU
                                # Velocity: km/s -> AU/day (1 day = 86400 seconds)
                                vx_au_per_day = (vx_kms * 86400) * KM_TO_AU
                                vy_au_per_day = (vy_kms * 86400) * KM_TO_AU
                                
                                # Verify these are reasonable AU values
                                if abs(x_au) < 100 and abs(y_au) < 100:  # Planets should be within 100 AU
                                    return (x_au, y_au, vx_au_per_day, vy_au_per_day)
                            except (ValueError, IndexError):
                                continue
        
        return None
    except (urllib.error.URLError, TimeoutError, ValueError, KeyError, IndexError) as e:
        # Silently fail and return None to trigger fallback
        return None


def fetch_planet_positions_online(epoch: dt.datetime) -> Optional[List[Tuple[str, str, float, float, float, float]]]:
    """
    Attempt to fetch real-time planet positions and velocities from NASA HORIZONS API.
    Returns list of (name, symbol, x, y, vx, vy) positions and velocities, or None if fetch fails.
    """
    positions: List[Tuple[str, str, float, float, float, float]] = []
    
    for name, _period_days, _radius, _lon0, symbol in PLANETS:
        coords = fetch_horizons_position(name, epoch)
        if coords is None:
            # If any planet fails, return None to trigger fallback
            return None
        x, y, vx, vy = coords
        positions.append((name, symbol, x, y, vx, vy))
    
    for name, symbol in ADDITIONAL_INNER_OBJECTS:
        coords = fetch_horizons_position(name, epoch)
        if coords is None:
            # Skip if unavailable; spacecraft ephemerides may be limited
            continue
        x, y, vx, vy = coords
        positions.append((name, symbol, x, y, vx, vy))
    
    for name, symbol in ADDITIONAL_OUTER_OBJECTS:
        coords = fetch_horizons_position(name, epoch)
        if coords is None:
            continue
        x, y, vx, vy = coords
        positions.append((name, symbol, x, y, vx, vy))
    
    # Add comet position using calculation (HORIZONS doesn't have 3I/ATLAS easily)
    # For comet, we'll estimate velocity from orbital mechanics or set to 0
    comet_position = compute_comet_position(epoch)
    if comet_position:
        name, symbol, x, y = comet_position
        # Estimate velocity for comet (simplified - could be improved)
        # For now, set velocities to 0 and let interpolation handle it
        positions.append((name, symbol, x, y, 0.0, 0.0))
    
    return positions


def interpolate_positions(
    base_positions: List[Tuple[str, str, float, float, float, float]],
    base_time: dt.datetime,
    current_time: dt.datetime,
) -> List[Tuple[str, str, float, float]]:
    """
    Interpolate positions based on velocity and elapsed time.
    Returns list of (name, symbol, x, y) interpolated positions.
    """
    elapsed_days = (current_time - base_time).total_seconds() / 86400.0
    interpolated = []
    
    for name, symbol, x_base, y_base, vx, vy in base_positions:
        # Linear interpolation: new_pos = base_pos + velocity * elapsed_time
        x_new = x_base + vx * elapsed_days
        y_new = y_base + vy * elapsed_days
        interpolated.append((name, symbol, x_new, y_new))
    
    return interpolated


def compute_planet_positions(epoch: dt.datetime, use_online: bool = True) -> List[Tuple[str, str, float, float]]:
    """
    Returns list of (name, symbol, x, y) positions in AU.
    Attempts to fetch real-time data from NASA HORIZONS API when online,
    falls back to calculated circular orbits if offline or fetch fails.
    """
    # Try online fetch first if enabled
    if use_online:
        online_positions = fetch_planet_positions_online(epoch)
        if online_positions is not None:
            return online_positions
    
    # Fallback to calculated positions using circular orbits
    days = _days_since_epoch(epoch)
    positions = []
    for name, period_days, radius, lon0, symbol in PLANETS:
        mean_motion = 360.0 / period_days
        angle_deg = (lon0 + mean_motion * days) % 360.0
        angle = math.radians(angle_deg)
        x = math.cos(angle) * radius
        y = math.sin(angle) * radius
        positions.append((name, symbol, x, y))
    comet_position = compute_comet_position(epoch)
    if comet_position:
        positions.append(comet_position)
    return positions


def scale_radius(
    radius: float,
    max_pixels: int,
    min_radius: float = MIN_ORBIT_RADIUS,
    max_radius: float = MAX_ORBIT_RADIUS,
) -> float:
    """
    Scale a planet's orbital radius to pixel coordinates using logarithmic spacing.

    The logarithmic mapping keeps inner planets visible while outer planets still
    occupy the edges of the viewport. Callers can override the minimum and maximum
    radius to re-use the helper for zoomed-in views (e.g., inner solar system).
    """

    if radius <= 0 or max_pixels <= 0:
        return 0.0

    min_radius = max(min_radius, 1e-6)
    max_radius = max(max_radius, min_radius)
    clamped_radius = max(radius, min_radius)

    log_min = math.log(min_radius)
    log_max = math.log(max_radius)
    log_radius = math.log(clamped_radius)

    if math.isclose(log_min, log_max):
        normalized_log = 1.0
    else:
        normalized_log = (log_radius - log_min) / (log_max - log_min)

    normalized = MIN_RADIUS_FRACTION + (1.0 - MIN_RADIUS_FRACTION) * normalized_log
    normalized = max(MIN_RADIUS_FRACTION, min(1.0, normalized))
    return normalized * max_pixels


def draw_orbit_path(win: curses.window, center_y: int, center_x: int, radius_px: float) -> None:
    """Render a dotted circle approximating an orbital path."""

    radius_int = int(round(radius_px))
    if radius_int <= 0:
        return

    steps = max(16, radius_int * 12)
    height, width = win.getmaxyx()
    for step in range(steps):
        angle = (2 * math.pi * step) / steps
        px = int(round(center_x + math.cos(angle) * radius_px))
        py = int(round(center_y - math.sin(angle) * radius_px * CHAR_ASPECT_RATIO))
        if 0 <= py < height and 0 <= px < width:
            try:
                win.addch(py, px, ".")
            except curses.error:
                continue


def render_planet_chart(
    win: curses.window,
    positions: List[Tuple[str, str, float, float]],
    title: str,
    min_radius: float,
    max_radius: float,
) -> None:
    """Render a scaled, color-aware planet map inside a curses subwindow."""

    win.erase()
    height, width = win.getmaxyx()
    if height < 4 or width < 6 or not positions:
        return

    # safe_addstr(win, 0, 0, title[: max(0, width - 1)])

    plot_top = 1
    plot_height = max(1, height - 2)
    center_y = plot_top + plot_height // 2
    center_x = width // 2
    vertical_limit = max(1, int((plot_height // 2) / CHAR_ASPECT_RATIO))
    max_radius_px = max(1, min(width // 2, vertical_limit))

    # Draw the sun marker at the chart centre and reserve those cells.
    occupied: Dict[Tuple[int, int], str] = {}
    sun_cells: set[Tuple[int, int]] = set()
    for dy, dx, char in [(0, 0, "O"), (0, -1, "-"), (0, 1, "-"), (-1, 0, "|"), (1, 0, "|")]:
        sy = center_y + dy
        sx = center_x + dx
        if 0 <= sy < height and 0 <= sx < width:
            try:
                win.addch(sy, sx, char)
            except curses.error:
                continue
            occupied[(sy, sx)] = "SUN"
            sun_cells.add((sy, sx))

    drawn_paths: set[int] = set()

    for name, symbol, x, y in positions:
        radius = math.hypot(x, y)
        radius_px = scale_radius(radius, max_radius_px, min_radius, max_radius)
        radius_int = int(round(radius_px))

        if radius_int > 1 and radius_int not in drawn_paths:
            draw_orbit_path(win, center_y, center_x, radius_px)
            drawn_paths.add(radius_int)

        angle = math.atan2(y, x)
        draw_x = int(round(center_x + math.cos(angle) * radius_px))
        draw_y = int(round(center_y - math.sin(angle) * radius_px * CHAR_ASPECT_RATIO))

        offsets = [(0, 0), (0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)]
        candidate_positions: List[Tuple[int, int]] = []
        seen_positions: set[Tuple[int, int]] = set()
        for distance in range(0, 3):  # search radius expands up to 2 cells away
            for dy, dx in offsets:
                candidate = (draw_y + dy * distance, draw_x + dx * distance)
                if candidate in seen_positions:
                    continue
                seen_positions.add(candidate)
                candidate_positions.append(candidate)
        
        placed = False
        for candidate_y, candidate_x in candidate_positions:
            if not (0 <= candidate_y < height and 0 <= candidate_x < width):
                continue
            if (candidate_y, candidate_x) in occupied:
                continue
            try:
                if COLOR_SUPPORT and name in PLANET_COLOR_PAIRS:
                    win.attron(curses.color_pair(PLANET_COLOR_PAIRS[name]) | curses.A_BOLD)
                    win.addch(candidate_y, candidate_x, symbol)
                    win.attroff(curses.color_pair(PLANET_COLOR_PAIRS[name]) | curses.A_BOLD)
                else:
                    win.addch(candidate_y, candidate_x, symbol)
            except curses.error:
                continue
            occupied[(candidate_y, candidate_x)] = name
            placed = True
            break

        if not placed:
            for candidate_y, candidate_x in candidate_positions:
                if not (0 <= candidate_y < height and 0 <= candidate_x < width):
                    continue
                if (candidate_y, candidate_x) in sun_cells:
                    continue
                try:
                    if COLOR_SUPPORT and name in PLANET_COLOR_PAIRS:
                        win.attron(curses.color_pair(PLANET_COLOR_PAIRS[name]) | curses.A_BOLD)
                        win.addch(candidate_y, candidate_x, symbol)
                        win.attroff(curses.color_pair(PLANET_COLOR_PAIRS[name]) | curses.A_BOLD)
                    else:
                        win.addch(candidate_y, candidate_x, symbol)
                except curses.error:
                    continue
                occupied[(candidate_y, candidate_x)] = name
                placed = True
                break

        if not placed:
            # As a last resort, drop the symbol at the computed position
            # this may overlap, but ensures the body is visible.
            base_y, base_x = draw_y, draw_x
            if 0 <= base_y < height and 0 <= base_x < width:
                try:
                    if COLOR_SUPPORT and name in PLANET_COLOR_PAIRS:
                        win.attron(curses.color_pair(PLANET_COLOR_PAIRS[name]) | curses.A_BOLD)
                        win.addch(base_y, base_x, symbol)
                        win.attroff(curses.color_pair(PLANET_COLOR_PAIRS[name]) | curses.A_BOLD)
                    else:
                        win.addch(base_y, base_x, symbol)
                    occupied[(base_y, base_x)] = name
                except curses.error:
                    pass


def compute_radius_bounds(
    positions: List[Tuple[str, str, float, float]]
) -> Tuple[float, float]:
    """Return min/max orbital radius (in AU) for a set of planet positions."""

    if not positions:
        return (MIN_ORBIT_RADIUS, MAX_ORBIT_RADIUS)

    radii = [max(1e-6, math.hypot(x, y)) for _, _, x, y in positions]
    min_radius = min(radii)
    max_radius = max(radii)
    if math.isclose(min_radius, max_radius):
        max_radius = max_radius * 1.1
    return (min_radius, max_radius)


def safe_addstr(win: curses.window, y: int, x: int, text: str) -> None:
    height, width = win.getmaxyx()
    if y >= height or x >= width:
        return
    trimmed = text[: max(0, width - x - 1)]
    if trimmed:
        try:
            win.addstr(y, x, trimmed)
        except curses.error:
            pass


def draw_window_box(win: curses.window, title: str) -> None:
    win.box()
    if len(title) > win.getmaxyx()[1] - 4:
        title = title[: max(0, win.getmaxyx()[1] - 4)]
    safe_addstr(win, 0, 2, f" {title} ")


def render_top_left(
    win: curses.window,
    weather_info: Optional[Dict[str, Any]],
    weather_error: Optional[str],
    animation_phase: int,
    power_bar_state: str,
) -> str:
    win.erase()
    header_text = "Local Weather"
    draw_window_box(win, header_text)

    height, width = win.getmaxyx()
    if height <= 4 or width <= 10:
        return power_bar_state

    content_width = max(1, width - 4)
    interior_height = max(0, height - 2)
    if interior_height <= 0:
        return power_bar_state

    power_bar_height = 0
    if interior_height >= POWER_BAR_MIN_HEIGHT + 2:
        desired_height = max(
            POWER_BAR_MIN_HEIGHT,
            min(POWER_BAR_MAX_HEIGHT, interior_height // 4),
        )
        desired_height = min(desired_height, interior_height - 3)
        if desired_height >= POWER_BAR_MIN_HEIGHT:
            power_bar_height = int(desired_height)

    weather_rows = max(0, interior_height - power_bar_height)
    weather_lines: List[str] = []

    if weather_info and weather_rows > 0:

        temp_c = weather_info.get("temperature_c")
        temp_f = weather_info.get("temperature_f")
        app_c = weather_info.get("apparent_c")
        app_f = weather_info.get("apparent_f")
        condition = weather_info.get("condition") or "--"

        wind_mph = weather_info.get("windspeed_mph")
        wind_kmh = weather_info.get("windspeed_kmh")
        wind_cardinal = weather_info.get("wind_direction_cardinal") or "--"

        humidity = weather_info.get("humidity")
        dew_c = weather_info.get("dewpoint_c")
        dew_f = weather_info.get("dewpoint_f")
        cloud_cover = weather_info.get("cloud_cover")
        visibility_km = weather_info.get("visibility_km")
        visibility_miles = weather_info.get("visibility_miles")
        surface_pressure_hpa = weather_info.get("surface_pressure_hpa")
        surface_pressure_inhg = weather_info.get("surface_pressure_inhg")
        pressure_msl_hpa = weather_info.get("pressure_msl_hpa")
        pressure_msl_inhg = weather_info.get("pressure_msl_inhg")
        precip_mm = weather_info.get("precipitation_mm")
        precip_in = weather_info.get("precipitation_in")
        rain_mm = weather_info.get("rain_mm")
        rain_in = weather_info.get("rain_in")
        snow_mm = weather_info.get("snow_mm")
        snow_in = weather_info.get("snow_in")
        uv_index_max = weather_info.get("uv_index_max")
        precip_probability = weather_info.get("precipitation_probability_max")
        sunrise_dt = weather_info.get("sunrise")
        sunset_dt = weather_info.get("sunset")
        is_day = weather_info.get("is_day")

        info_lines: List[str] = [f"Conditions: {condition}"]

        if isinstance(dew_c, (int, float)) and isinstance(dew_f, (int, float)):
            info_lines.append(f"Dew point: {dew_c:.1f}°C ({dew_f:.1f}°F)")

        if isinstance(cloud_cover, (int, float)):
            info_lines.append(f"Clouds: {cloud_cover:.0f}% cover")

        if isinstance(visibility_km, (int, float)):
            if isinstance(visibility_miles, (int, float)):
                info_lines.append(
                    f"Visibility: {visibility_km:.1f} km ({visibility_miles:.1f} mi)"
                )
            else:
                info_lines.append(f"Visibility: {visibility_km:.1f} km")

        if isinstance(surface_pressure_hpa, (int, float)):
            pressure_line = f"Pressure: {surface_pressure_hpa:.0f} hPa"
            if isinstance(surface_pressure_inhg, (int, float)):
                pressure_line += f" ({surface_pressure_inhg:.2f} inHg)"
            info_lines.append(pressure_line)

        if isinstance(pressure_msl_hpa, (int, float)):
            msl_line = f"MSL Pressure: {pressure_msl_hpa:.0f} hPa"
            if isinstance(pressure_msl_inhg, (int, float)):
                msl_line += f" ({pressure_msl_inhg:.2f} inHg)"
            info_lines.append(msl_line)

        def format_precip_line(label: str, mm_value: Optional[float], inch_value: Optional[float]) -> None:
            if mm_value is None:
                return
            line = f"{label}: {mm_value:.1f} mm"
            if inch_value is not None:
                line += f" ({inch_value:.2f} in)"
            info_lines.append(line)

        format_precip_line("Precip (1h)", precip_mm, precip_in)

        if isinstance(rain_mm, (int, float)) and rain_mm > 0:
            format_precip_line(
                "Rain (1h)",
                rain_mm,
                rain_in if isinstance(rain_in, (int, float)) else None,
            )

        if isinstance(snow_mm, (int, float)) and snow_mm > 0:
            format_precip_line(
                "Snow (1h)",
                snow_mm,
                snow_in if isinstance(snow_in, (int, float)) else None,
            )

        if isinstance(sunrise_dt, dt.datetime):
            info_lines.append(f"Sunrise: {sunrise_dt.strftime('%H:%M UTC')}")

        if isinstance(sunset_dt, dt.datetime):
            info_lines.append(f"Sunset: {sunset_dt.strftime('%H:%M UTC')}")

        if isinstance(is_day, bool):
            info_lines.append(f"Daylight: {'Yes' if is_day else 'No'}")

        gauge_lines: List[str] = []

        if isinstance(temp_c, (int, float)):
            temp_display = (
                f"{temp_c:.1f}°C / {temp_f:.1f}°F"
                if isinstance(temp_f, (int, float))
                else f"{temp_c:.1f}°C"
            )
            gauge_line = make_weather_gauge_line(
                "Temp",
                float(temp_c),
                -30.0,
                45.0,
                temp_display,
                content_width,
                animation_phase,
            )
            if gauge_line:
                gauge_lines.append(gauge_line)

        if isinstance(app_c, (int, float)):
            feels_display = (
                f"{app_c:.1f}°C / {app_f:.1f}°F"
                if isinstance(app_f, (int, float))
                else f"{app_c:.1f}°C"
            )
            gauge_line = make_weather_gauge_line(
                "Feels",
                float(app_c),
                -30.0,
                45.0,
                feels_display,
                content_width,
                animation_phase,
            )
            if gauge_line:
                gauge_lines.append(gauge_line)

        wind_gauge_value: Optional[float] = None
        if isinstance(wind_mph, (int, float)):
            wind_gauge_value = float(wind_mph)
        elif isinstance(wind_kmh, (int, float)):
            wind_gauge_value = float(wind_kmh) * 0.621371
        if wind_gauge_value is not None:
            wind_parts: List[str] = []
            if isinstance(wind_mph, (int, float)):
                wind_parts.append(f"{wind_mph:.1f} mph")
            if isinstance(wind_kmh, (int, float)):
                wind_parts.append(f"{wind_kmh:.1f} km/h")
            wind_display = " / ".join(wind_parts) or f"{wind_gauge_value:.1f} mph"
            if wind_cardinal and wind_cardinal != "--":
                wind_display = f"{wind_display} {wind_cardinal}"
            gauge_line = make_weather_gauge_line(
                "Wind",
                wind_gauge_value,
                0.0,
                60.0,
                wind_display,
                content_width,
                animation_phase,
            )
            if gauge_line:
                gauge_lines.append(gauge_line)

        if isinstance(humidity, (int, float)):
            gauge_line = make_weather_gauge_line(
                "Humidity",
                float(humidity),
                0.0,
                100.0,
                f"{humidity:.0f}%",
                content_width,
                animation_phase,
            )
            if gauge_line:
                gauge_lines.append(gauge_line)

        if isinstance(precip_probability, (int, float)):
            gauge_line = make_weather_gauge_line(
                "Precip prob",
                float(precip_probability),
                0.0,
                100.0,
                f"{precip_probability:.0f}%",
                content_width,
                animation_phase,
            )
            if gauge_line:
                gauge_lines.append(gauge_line)

        if isinstance(uv_index_max, (int, float)):
            gauge_line = make_weather_gauge_line(
                "UV",
                float(uv_index_max),
                0.0,
                11.0,
                f"{uv_index_max:.1f}",
                content_width,
                animation_phase,
            )
            if gauge_line:
                gauge_lines.append(gauge_line)

        if weather_lines and (gauge_lines or info_lines):
            weather_lines.append("")

        if gauge_lines:
            weather_lines.extend(gauge_lines)

        if gauge_lines and info_lines:
            weather_lines.append("")

        weather_lines.extend(info_lines)

        if weather_error:
            weather_lines.append("")
            weather_lines.extend(wrap_text_lines(weather_error, content_width))

    else:
        message = weather_error or (
            f"Set {WEATHER_LAT_ENV}/{WEATHER_LON_ENV} env vars to enable weather."
        )
        if weather_rows > 0:
            weather_lines.extend(wrap_text_lines(message, content_width))

    while weather_lines and weather_lines[0] == "":
        weather_lines.pop(0)
    while weather_lines and weather_lines[-1] == "":
        weather_lines.pop()

    row = 1
    rows_drawn = 0
    for line in weather_lines:
        if weather_rows and rows_drawn >= weather_rows:
            break
        if row >= height - 1:
            break
        if line:
            safe_addstr(win, row, 2, line[: content_width])
        row += 1
        rows_drawn += 1

    if power_bar_height > 0:
        separator_row = min(
            height - power_bar_height - 1,
            max(1, 1 + (weather_rows if weather_rows > 0 else rows_drawn)),
        )
        if separator_row < height - 1:
            for x in range(1, width - 1):
                try:
                    win.addch(separator_row, x, "-")
                except curses.error:
                    break
            try:
                win.addch(separator_row, 0, "+")
                win.addch(separator_row, width - 1, "+")
            except curses.error:
                pass

        power_start_row = min(height - power_bar_height, separator_row + 1)
        max_rows_possible = max(0, height - power_start_row - 1)
        available_rows = min(power_bar_height, max_rows_possible)
        if available_rows > 0:
            power_bar_state = render_power_bar_section(
                win,
                power_start_row,
                available_rows,
                content_width,
                power_bar_state,
                animation_phase,
            )

    return power_bar_state


def render_power_bar_section(
    win: curses.window,
    start_row: int,
    available_rows: int,
    content_width: int,
    buffer: str,
    animation_phase: int,
) -> str:
    if available_rows <= 0 or content_width <= 0:
        return buffer

    height, _ = win.getmaxyx()
    current_row = start_row
    rows_left = available_rows

    spinner = WEATHER_SPINNER_CHARS[animation_phase % len(WEATHER_SPINNER_CHARS)]
    status_line = f"Power Feed {spinner}"
    if current_row < height - 1:
        safe_addstr(win, current_row, 2, status_line[: content_width])
    current_row += 1
    rows_left -= 1

    if rows_left <= 0 or current_row >= height - 1:
        return buffer

    bar_width = max(8, content_width)
    buffer = update_power_bar_buffer(buffer, bar_width)
    bar_text = f"[{buffer}]"
    safe_addstr(win, current_row, 2, bar_text[: content_width])
    current_row += 1
    rows_left -= 1

    if rows_left <= 0 or current_row >= height - 1:
        return buffer

    active_cells = sum(1 for ch in buffer if ch not in (" ", "."))
    utilization = int((active_cells / len(buffer)) * 100) if buffer else 0
    output_line = f"Output: {utilization:3d}%"
    safe_addstr(win, current_row, 2, output_line[: content_width])

    return buffer


def render_bottom_left(
    win: curses.window,
    history: List[str],
    partial_line: str,
    last_refresh: dt.datetime | None,
) -> None:
    win.erase()
    header_text = "Live News Feeds (Breaking News)"

    draw_window_box(win, header_text)

    height, width = win.getmaxyx()
    if height <= 5 or width <= 10:
        return

    content_width = max(1, width - 4)
    feed_rows = max(0, height - 2)

    if feed_rows <= 0:
        return

    wrapped_lines: List[str] = []
    for line in history:
        wrapped_lines.extend(wrap_text_lines(line, content_width))

    if partial_line:
        wrapped_lines.extend(wrap_text_lines(partial_line, content_width))

    if not wrapped_lines:
        placeholder = "Awaiting latest headlines..."
        safe_addstr(win, 1, 2, placeholder[: max(0, width - 4)])
        return

    visible_lines = wrapped_lines[-feed_rows:]
    start_row = max(1, height - 1 - len(visible_lines))

    for offset, text in enumerate(visible_lines):
        y = start_row + offset
        if y >= height - 1:
            break
        safe_addstr(win, y, 2, text[: max(0, width - 4)])


def render_telemetry_panel(
    win: curses.window,
    positions: List[Tuple[str, str, float, float]],
    title: str = "Telemetry",
    previous_telemetry: Optional[Dict[str, Tuple[float, float, float, float]]] = None,
    change_times: Optional[Dict[str, float]] = None,
    current_time: Optional[float] = None,
    highlight_duration: float = 1.0,
) -> None:
    """
    Render a telemetry panel showing position and distance from Sun for each object.
    Highlights changed values in bold for a brief period.
    """
    win.erase()
    draw_window_box(win, title)
    
    height, width = win.getmaxyx()
    if height <= 3 or width <= 10:
        return
    
    # Sort positions by distance from sun (closest first)
    positions_with_distance = []
    for name, symbol, x, y in positions:
        distance = math.hypot(x, y)
        positions_with_distance.append((name, symbol, x, y, distance))
    positions_with_distance.sort(key=lambda p: p[4])  # Sort by distance
    
    content_width = max(1, width - 4)
    start_row = 1
    max_rows = height - 2
    
    # Track changes: previous_telemetry stores (prev_x, prev_y, prev_dist, last_update_time)
    # change_times tracks when each individual value (x, y, dist) last changed
    # Update change_times dict in-place when values change
    if previous_telemetry is not None and current_time is not None and change_times is not None:
        for name, symbol, x, y, distance in positions_with_distance:
            key = name
            if key in previous_telemetry:
                prev_x, prev_y, prev_dist, _ = previous_telemetry[key]
                # Check if values changed (with small tolerance for floating point)
                x_changed = not math.isclose(x, prev_x, abs_tol=1e-8)
                y_changed = not math.isclose(y, prev_y, abs_tol=1e-8)
                dist_changed = not math.isclose(distance, prev_dist, abs_tol=1e-8)
                
                # Track when each value changed - update the passed-in dict
                if x_changed:
                    change_times[f"{key}_x"] = current_time
                if y_changed:
                    change_times[f"{key}_y"] = current_time
                if dist_changed:
                    change_times[f"{key}_dist"] = current_time
                
                # Update stored values for next comparison
                previous_telemetry[key] = (x, y, distance, current_time)
            else:
                # First time seeing this object - store initial values
                previous_telemetry[key] = (x, y, distance, current_time)
                # Mark all as just changed (first appearance)
                change_times[f"{key}_x"] = current_time
                change_times[f"{key}_y"] = current_time
                change_times[f"{key}_dist"] = current_time
    
    for idx, (name, symbol, x, y, distance) in enumerate(positions_with_distance):
        if idx >= max_rows:
            break
        
        row = start_row + idx
        if row >= height - 1:
            break
        
        # Format: symbol name: X=... Y=... Dist=...AU
        # Use higher precision to show small changes (6 decimal places)
        x_str = f"{x:+.6f}"
        y_str = f"{y:+.6f}"
        dist_str = f"{distance:.6f}"
        
        # Check if values should be highlighted (recently changed)
        x_highlight = False
        y_highlight = False
        dist_highlight = False
        
        if current_time is not None and change_times is not None:
            key = name
            x_key = f"{key}_x"
            y_key = f"{key}_y"
            dist_key = f"{key}_dist"
            
            if x_key in change_times:
                time_since_change = current_time - change_times[x_key]
                x_highlight = time_since_change < highlight_duration and time_since_change >= 0
            
            if y_key in change_times:
                time_since_change = current_time - change_times[y_key]
                y_highlight = time_since_change < highlight_duration and time_since_change >= 0
            
            if dist_key in change_times:
                time_since_change = current_time - change_times[dist_key]
                dist_highlight = time_since_change < highlight_duration and time_since_change >= 0
        
        # Build telemetry line with conditional bold formatting
        prefix = f"{symbol} {name[:8]:8s} "
        x_part = f"X:{x_str:>12s} "
        y_part = f"Y:{y_str:>12s} "
        dist_part = f"D:{dist_str:>10s}AU"
        
        # Write prefix
        col = 2
        safe_addstr(win, row, col, prefix[:content_width])
        col += len(prefix)
        
        # Write X with optional bold
        if x_highlight and COLOR_SUPPORT:
            try:
                win.attron(curses.A_BOLD)
                win.addstr(row, col, x_part[:content_width - (col - 2)])
                win.attroff(curses.A_BOLD)
            except curses.error:
                safe_addstr(win, row, col, x_part[:content_width - (col - 2)])
        else:
            safe_addstr(win, row, col, x_part[:content_width - (col - 2)])
        col += len(x_part)
        
        # Write Y with optional bold
        if y_highlight and COLOR_SUPPORT:
            try:
                win.attron(curses.A_BOLD)
                win.addstr(row, col, y_part[:content_width - (col - 2)])
                win.attroff(curses.A_BOLD)
            except curses.error:
                safe_addstr(win, row, col, y_part[:content_width - (col - 2)])
        else:
            safe_addstr(win, row, col, y_part[:content_width - (col - 2)])
        col += len(y_part)
        
        # Write distance with optional bold
        if dist_highlight and COLOR_SUPPORT:
            try:
                win.attron(curses.A_BOLD)
                win.addstr(row, col, dist_part[:content_width - (col - 2)])
                win.attroff(curses.A_BOLD)
            except curses.error:
                safe_addstr(win, row, col, dist_part[:content_width - (col - 2)])
        else:
            safe_addstr(win, row, col, dist_part[:content_width - (col - 2)])


def render_bottom_right(
    win: curses.window,
    positions: List[Tuple[str, str, float, float]],
) -> None:
    win.erase()
    draw_window_box(win, "Outer System (top-down) - Gas & ice giants perspective")

    height, width = win.getmaxyx()
    if height <= 4 or width <= 10:
        return

    outer_positions = [pos for pos in positions if pos[0] in OUTER_PLANET_NAMES]
    if not outer_positions:
        safe_addstr(win, height // 2, 2, "Outer-planet data unavailable")
        return

    chart_height = max(4, height - 2)
    chart_width = max(4, width - 2)
    chart_window = win.derwin(chart_height, chart_width, 1, 1)
    outer_min, outer_max = compute_radius_bounds(outer_positions)
    render_planet_chart(
        chart_window,
        outer_positions,
        "Outer system (Jupiter–Neptune)",
        outer_min * 0.9,
        outer_max * 1.1,
    )


def render_top_right(
    win: curses.window,
    positions: List[Tuple[str, str, float, float]],
    last_refresh: dt.datetime | None,
    seconds_since_refresh: float,
    data_source: str = "calculated",
) -> None:
    win.erase()
    refresh_text = (
        last_refresh.strftime("Last refresh %H:%M:%S UTC") if last_refresh else "Waiting for first update..."
    )
    next_update_text = f"Next sync update in {int(max(0, SYNC_UPDATE_INTERVAL_SECONDS - seconds_since_refresh))}s"
    source_text = f"[{data_source}]"
    header_text = f"Solar System (top-down) - {refresh_text} - {next_update_text} - {source_text}"
    draw_window_box(win, header_text)

    height, width = win.getmaxyx()
    if height <= 4 or width <= 10:
        return

    inner_positions = [pos for pos in positions if pos[0] in INNER_PLANET_NAMES]
    if not inner_positions:
        safe_addstr(win, height // 2, 2, "Inner-planet data unavailable")
        return

    chart_height = max(4, height - 2)
    chart_width = max(4, width - 2)
    inner_window = win.derwin(chart_height, chart_width, 1, 1)
    inner_min, inner_max = compute_radius_bounds(inner_positions)
    render_planet_chart(
        inner_window,
        inner_positions,
        "Inner system (Mercury–Mars)",
        inner_min * 0.95,
        inner_max * 1.25,
    )


def draw_dashboard(
    stdscr: curses.window,
    positions: List[Tuple[str, str, float, float]],
    last_refresh: dt.datetime | None,
    seconds_since_refresh: float,
    weather_info: Optional[Dict[str, Any]],
    weather_error: Optional[str],
    weather_animation_phase: int,
    power_bar_state: str,
    news_history: List[str],
    news_partial_line: str,
    headlines_updated: dt.datetime | None,
    data_source: str = "calculated",
    previous_telemetry: Optional[Dict[str, Tuple[float, float, float, float]]] = None,
    change_times: Optional[Dict[str, float]] = None,
    current_time: Optional[float] = None,
    highlight_duration: float = 1.0,
) -> str:
    stdscr.erase()
    height, width = stdscr.getmaxyx()
    min_height, min_width = 16, 60
    if height < min_height or width < min_width:
        warning = "Terminal too small. Resize to at least 60x16."
        stdscr.addstr(0, 0, warning[: width - 1])
        stdscr.refresh()
        return power_bar_state

    left_width = max(20, int(width * 0.25))
    telemetry_width = max(20, int(width * 0.25))
    mid_width = width - left_width - telemetry_width
    if mid_width < 20:
        # ensure the central chart column remains usable
        deficit = 20 - mid_width
        if telemetry_width - deficit >= 20:
            telemetry_width -= deficit
        else:
            remaining_deficit = deficit - (telemetry_width - 20)
            telemetry_width = 20
            left_width = max(20, left_width - remaining_deficit)
        mid_width = width - left_width - telemetry_width

    top_height = height // 2
    bottom_height = height - top_height

    windows = {
        "weather": stdscr.derwin(top_height, left_width, 0, 0),
        "news": stdscr.derwin(bottom_height, left_width, top_height, 0),
        "inner_chart": stdscr.derwin(top_height, mid_width, 0, left_width),
        "outer_chart": stdscr.derwin(bottom_height, mid_width, top_height, left_width),
        "inner_tel": stdscr.derwin(top_height, telemetry_width, 0, left_width + mid_width),
        "outer_tel": stdscr.derwin(bottom_height, telemetry_width, top_height, left_width + mid_width),
    }

    power_bar_state = render_top_left(
        windows["weather"],
        weather_info,
        weather_error,
        weather_animation_phase,
        power_bar_state,
    )
    render_top_right(
        windows["inner_chart"],
        positions,
        last_refresh,
        seconds_since_refresh,
        data_source,
    )
    render_bottom_left(windows["news"], news_history, news_partial_line, headlines_updated)
    render_bottom_right(windows["outer_chart"], positions)

    inner_positions = [pos for pos in positions if pos[0] in INNER_PLANET_NAMES]
    outer_positions = [pos for pos in positions if pos[0] in OUTER_PLANET_NAMES]
    if positions:
        extras = [pos for pos in positions if pos[0] not in INNER_PLANET_NAMES and pos[0] not in OUTER_PLANET_NAMES]
        if extras:
            outer_positions = outer_positions + extras

    render_telemetry_panel(
        windows["inner_tel"],
        inner_positions,
        "Inner System Telemetry",
        previous_telemetry,
        change_times,
        current_time,
        highlight_duration,
    )
    render_telemetry_panel(
        windows["outer_tel"],
        outer_positions,
        "Outer System Telemetry",
        previous_telemetry,
        change_times,
        current_time,
        highlight_duration,
    )
 
    stdscr.refresh()
    return power_bar_state


def main(stdscr: curses.window) -> None:
    curses.curs_set(0)
    init_color_support()
    load_env_files() # Load environment variables from .env files
    stdscr.nodelay(True)
    stdscr.timeout(250)

    # Synchronized update cycle: all components update together every 5 minutes
    last_sync_update_ts = 0.0
    last_refresh_dt: dt.datetime | None = None
    positions: List[Tuple[str, str, float, float]] = []
    # Store base positions with velocities for interpolation
    base_positions_with_velocities: List[Tuple[str, str, float, float, float, float]] = []
    base_positions_timestamp: dt.datetime | None = None
    last_telemetry_update_ts = 0.0
    TELEMETRY_UPDATE_INTERVAL = 0.5  # Update telemetry every 0.5 seconds

    news_items: List[str] = []
    news_pointer = 0
    news_last_update_dt: dt.datetime | None = None
    news_history: List[str] = []
    news_current: str = ""
    news_progress = 0
    news_char_ts = 0.0

    weather_coords: Optional[Tuple[float, float]] = load_weather_coordinates()
    weather_data: Optional[Dict[str, Any]] = None
    weather_error: Optional[str] = None
    weather_animation_phase = 0
    power_bar_state = ""
    data_source = "calculated"  # Track whether using online or calculated data
    
    # Track previous telemetry values for change detection and highlighting
    previous_telemetry: Dict[str, Tuple[float, float, float, float]] = {}  # name -> (x, y, dist, timestamp)
    # Track when each value last changed: "name_x" -> timestamp, "name_y" -> timestamp, etc.
    telemetry_change_times: Dict[str, float] = {}
    TELEMETRY_HIGHLIGHT_DURATION = 1.0  # Highlight changed values for 1 second

    while True:
        now = dt.datetime.now(dt.timezone.utc)
        now_ts = now.timestamp()
        time_since_sync = now_ts - last_sync_update_ts
        
        # Synchronized update: all components refresh together every 5 minutes
        if not positions or time_since_sync >= SYNC_UPDATE_INTERVAL_SECONDS:
            last_refresh_dt = now
            last_sync_update_ts = now_ts
            
            # Update planet positions (try online first, fallback to calculated)
            try:
                online_positions = fetch_planet_positions_online(last_refresh_dt)
                if online_positions is not None:
                    # Store positions with velocities for interpolation
                    base_positions_with_velocities = online_positions
                    base_positions_timestamp = last_refresh_dt
                    # Initial positions (will be interpolated on next render)
                    positions = interpolate_positions(online_positions, last_refresh_dt, now)
                    data_source = "NASA HORIZONS"
                else:
                    # Fallback: no velocities available, use calculated positions
                    positions = compute_planet_positions(last_refresh_dt, use_online=False)
                    base_positions_with_velocities = []
                    base_positions_timestamp = None
                    data_source = "calculated"
            except Exception:
                # If anything fails, use calculated positions
                positions = compute_planet_positions(last_refresh_dt, use_online=False)
                base_positions_with_velocities = []
                base_positions_timestamp = None
                data_source = "calculated"
            
            last_telemetry_update_ts = now_ts
            
            # Update weather data
            current_coords = load_weather_coordinates()
            if current_coords is None:
                weather_coords = None
                weather_data = None
                weather_error = (
                    f"Configure {WEATHER_LAT_ENV}/{WEATHER_LON_ENV} to show local weather."
                )
            else:
                weather_error = None
                if current_coords != weather_coords:
                    weather_coords = current_coords
                lat, lon = current_coords
                fetched_weather = fetch_weather_data(lat, lon)
                if fetched_weather:
                    weather_data = fetched_weather
                    weather_error = None
                else:
                    weather_error = "Weather lookup failed (network or API issue)."
            
            # Update news headlines
            fetched = fetch_latest_headlines()
            if fetched:
                news_items = fetched
                news_pointer = 0
                news_current = ""
                news_progress = 0
                news_last_update_dt = now
        
        elapsed = time_since_sync

        # Interpolate positions continuously if we have velocity data
        # This makes the telemetry update smoothly in real-time
        if base_positions_with_velocities and base_positions_timestamp:
            # Always interpolate to current time for smooth updates
            positions = interpolate_positions(base_positions_with_velocities, base_positions_timestamp, now)
        elif not positions:
            # If we don't have positions at all, try to get calculated ones
            if last_refresh_dt:
                positions = compute_planet_positions(last_refresh_dt, use_online=False)

        # News ticker animation (independent of sync cycle)
        if news_items and not news_current:
            news_current = news_items[news_pointer % len(news_items)]
            news_pointer += 1
            news_progress = 0
            news_char_ts = now_ts

        if news_current and now_ts - news_char_ts >= NEWS_CHAR_INTERVAL:
            elapsed_chars = max(1, int((now_ts - news_char_ts) / NEWS_CHAR_INTERVAL))
            news_progress = min(len(news_current), news_progress + elapsed_chars)
            news_char_ts = now_ts
            if news_progress >= len(news_current):
                news_history.append(news_current)
                if len(news_history) > NEWS_HISTORY_LIMIT:
                    news_history = news_history[-NEWS_HISTORY_LIMIT:]
                news_current = ""
                news_progress = 0
                news_char_ts = now_ts

        weather_animation_phase = (weather_animation_phase + 1) % WEATHER_ANIMATION_MODULO

        power_bar_state = draw_dashboard(
            stdscr,
            positions,
            last_refresh_dt,
            elapsed,
            weather_data,
            weather_error,
            weather_animation_phase,
            power_bar_state,
            news_history,
            (news_current[:news_progress] + "_") if news_current else "",
            news_last_update_dt,
            data_source,
            previous_telemetry,
            telemetry_change_times,
            now_ts,
            TELEMETRY_HIGHLIGHT_DURATION,
        )

        key = stdscr.getch()
        if key in (ord("q"), ord("Q")):
            break
        if key in (ord("r"), ord("R")):
            last_sync_update_ts = 0.0  # trigger immediate sync update next loop


if __name__ == "__main__":
    curses.wrapper(main)
