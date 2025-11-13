/**
 * Space Operations Alpha Assistant Dashboard - JavaScript Implementation
 * 
 * This module provides a web-based dashboard that mirrors the functionality
 * of the Python curses-based terminal dashboard. It displays real-time solar
 * system positions, weather data, news feeds, and telemetry in a four-panel
 * layout. The dashboard fetches data from NASA HORIZONS API, Open-Meteo API,
 * and RSS news feeds, with smooth animations and interpolated position updates.
 */

// Configuration constants
const SYNC_UPDATE_INTERVAL_SECONDS = 300; // 5 minutes
const TELEMETRY_UPDATE_INTERVAL = 0.5; // 0.5 seconds
const TELEMETRY_HIGHLIGHT_DURATION = 1.0; // 1 second
const NEWS_CHAR_INTERVAL = 0.05; // 50ms per character
const NEWS_HISTORY_LIMIT = 120;
const NEWS_MAX_ITEMS = 40;
const CHAR_ASPECT_RATIO = 0.5;
const MIN_RADIUS_FRACTION = 0.05;
const J2000 = new Date('2000-01-01T12:00:00Z');

// Planet data: (name, period_days, radius_AU, mean_longitude_deg, symbol)
const PLANETS = [
    ["Mercury", 87.9691, 0.39, 252.250906, "m"],
    ["Venus", 224.701, 0.72, 181.979801, "V"],
    ["Earth", 365.256, 1.00, 100.466457, "E"],
    ["Mars", 686.980, 1.52, 355.433000, "M"],
    ["Jupiter", 4332.589, 5.20, 34.351519, "J"],
    ["Saturn", 10759.22, 9.58, 50.077444, "S"],
    ["Uranus", 30685.4, 19.20, 314.055005, "U"],
    ["Neptune", 60189.0, 30.05, 304.348665, "N"],
];

const ADDITIONAL_INNER_OBJECTS = [
    ["Parker Solar Probe", "p"],
    ["Juice Spacecraft", "j"],
];

const ADDITIONAL_OUTER_OBJECTS = [
    ["Europa Clipper", "c"],
];

const INNER_PLANET_NAMES = new Set([
    "Mercury", "Venus", "Earth", "Mars", "3I/ATLAS",
    "Parker Solar Probe", "Juice Spacecraft"
]);

const OUTER_PLANET_NAMES = new Set([
    "Jupiter", "Saturn", "Uranus", "Neptune", "3I/ATLAS",
    "Europa Clipper"
]);

const HORIZONS_PLANET_IDS = {
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
};

const HORIZONS_API_URL = "https://ssd.jpl.nasa.gov/api/horizons.api";
const HORIZONS_TIMEOUT = 10000;
// Use CORS proxy for browser requests (NASA API doesn't allow direct browser access)
const USE_CORS_PROXY = true;
const CORS_PROXY_URL = "https://api.allorigins.win/get?url=";

const WEATHER_API_URL = "https://api.open-meteo.com/v1/forecast";
const WEATHER_LAT_ENV = "ALPHA_WEATHER_LAT";
const WEATHER_LON_ENV = "ALPHA_WEATHER_LON";

const WEATHER_CODE_MAP = {
    0: "Clear sky", 1: "Mainly clear", 2: "Partly cloudy", 3: "Overcast",
    45: "Fog", 48: "Depositing rime fog", 51: "Light drizzle", 53: "Moderate drizzle",
    55: "Dense drizzle", 56: "Light freezing drizzle", 57: "Dense freezing drizzle",
    61: "Slight rain", 63: "Moderate rain", 65: "Heavy rain",
    66: "Light freezing rain", 67: "Heavy freezing rain",
    71: "Slight snow", 73: "Moderate snow", 75: "Heavy snow", 77: "Snow grains",
    80: "Light rain showers", 81: "Moderate rain showers", 82: "Violent rain showers",
    85: "Light snow showers", 86: "Heavy snow showers",
    95: "Thunderstorm", 96: "Thunderstorm with slight hail", 99: "Thunderstorm with heavy hail",
};

const NEWS_FEEDS = [
    ["NASA", "https://www.nasa.gov/rss/dyn/breaking_news.rss"],
    ["ESA", "https://www.esa.int/rssfeed/Our_Activities/Space_News"],
    ["SpaceNews", "https://spacenews.com/feed/"],
];

const PLANET_COLORS = {
    "Mercury": "#ffffff",
    "Venus": "#ffff00",
    "Earth": "#00ff00",
    "Mars": "#ff0000",
    "Jupiter": "#ff00ff",
    "Saturn": "#00ffff",
    "Uranus": "#0000ff",
    "Neptune": "#ffffff",
    "3I/ATLAS": "#ffff00",
    "Parker Solar Probe": "#ffff00",
    "Juice Spacecraft": "#00ff00",
    "Europa Clipper": "#00ffff",
};

const WEATHER_SPINNER_CHARS = ["-", "\\", "|", "/"];
const POWER_BAR_TEXTURES = [".", ":", "-", "=", "+", "*", "#", "%", "@"];

const GAUSSIAN_GRAVITATIONAL_CONSTANT = 0.01720209895;

const COMET_3I_ATLAS = [
    "3I/ATLAS",
    -0.2639,  // semi-major axis (AU, hyperbolic => negative)
    6.1396,   // eccentricity
    175.1131, // inclination (deg)
    322.157,  // longitude of ascending node (deg)
    128.010,  // argument of perihelion (deg)
    2460977.98, // time of perihelion passage (Julian Date)
    1.3564,   // perihelion distance (AU)
    "*",      // display symbol
];

// Global state
let lastSyncUpdateTs = 0;
let lastRefreshDt = null;
let positions = [];
let basePositionsWithVelocities = [];
let basePositionsTimestamp = null;
let dataSource = "calculated";

let weatherData = null;
let weatherError = null;
let weatherAnimationPhase = 0;
let powerBarState = "";
let lastPowerBarUpdate = 0;
const POWER_BAR_UPDATE_INTERVAL = 200; // Update power bar every 200ms

let newsItems = [];
let newsHistory = [];
let newsCurrent = "";
let newsProgress = 0;
let newsCharTs = 0;
let newsLastUpdateDt = null;

let previousTelemetry = {}; // name -> {x, y, dist, timestamp}
let telemetryChangeTimes = {}; // "name_x" -> timestamp, etc.

// Canvas contexts
let innerChartCtx = null;
let outerChartCtx = null;

/**
 * Initialize the dashboard on page load
 */
function initDashboard() {
    // Get canvas contexts
    const innerCanvas = document.getElementById('inner-chart-canvas');
    const outerCanvas = document.getElementById('outer-chart-canvas');
    
    if (innerCanvas) {
        innerChartCtx = innerCanvas.getContext('2d');
        resizeCanvas(innerCanvas);
    }
    
    if (outerCanvas) {
        outerChartCtx = outerCanvas.getContext('2d');
        resizeCanvas(outerCanvas);
    }
    
    // Handle window resize
    window.addEventListener('resize', () => {
        resizeCanvas(innerCanvas);
        resizeCanvas(outerCanvas);
    });
    
    // Start main loop
    requestAnimationFrame(mainLoop);
}

/**
 * Resize canvas to match container size
 */
function resizeCanvas(canvas) {
    const container = canvas.parentElement;
    if (container) {
        const rect = container.getBoundingClientRect();
        canvas.width = rect.width;
        canvas.height = rect.height;
    }
}

/**
 * Convert datetime to Julian Date
 */
function datetimeToJulianDate(dt) {
    const daysSinceJ2000 = (dt - J2000) / (1000 * 60 * 60 * 24);
    return 2451545.0 + daysSinceJ2000;
}

/**
 * Convert datetime to days since J2000
 */
function daysSinceEpoch(dt) {
    return (dt - J2000) / (1000 * 60 * 60 * 24);
}

/**
 * Solve hyperbolic anomaly using Newton-Raphson
 */
function solveHyperbolicAnomaly(meanAnomaly, eccentricity) {
    if (Math.abs(meanAnomaly) < 1e-10) return 0.0;
    if (eccentricity <= 1.0) throw new Error("Hyperbolic anomaly requires e > 1");
    
    let F = Math.asinh(meanAnomaly / eccentricity);
    for (let i = 0; i < 60; i++) {
        const sinhF = Math.sinh(F);
        const coshF = Math.cosh(F);
        const fVal = eccentricity * sinhF - F - meanAnomaly;
        const fPrime = eccentricity * coshF - 1.0;
        if (Math.abs(fPrime) < 1e-12) break;
        const delta = fVal / fPrime;
        F -= delta;
        if (Math.abs(delta) < 1e-10) break;
    }
    return F;
}

/**
 * Compute comet 3I/ATLAS position
 */
function computeCometPosition(epoch) {
    const [name, semiMajorAxis, eccentricity, inclinationDeg, longAscNodeDeg,
           argPerihelionDeg, perihelionJd, , symbol] = COMET_3I_ATLAS;
    
    const jd = datetimeToJulianDate(epoch);
    const deltaT = jd - perihelionJd;
    
    if (semiMajorAxis >= 0 || eccentricity <= 1) return null;
    
    const meanMotion = GAUSSIAN_GRAVITATIONAL_CONSTANT / Math.sqrt(Math.abs(semiMajorAxis) ** 3);
    const meanAnomaly = meanMotion * deltaT;
    
    let hyperbolicAnomaly;
    try {
        hyperbolicAnomaly = solveHyperbolicAnomaly(meanAnomaly, eccentricity);
    } catch (e) {
        return null;
    }
    
    const coshF = Math.cosh(hyperbolicAnomaly);
    const sinhF = Math.sinh(hyperbolicAnomaly);
    
    let radius = semiMajorAxis * (eccentricity * coshF - 1.0);
    radius = Math.abs(radius);
    
    const trueAnomaly = 2.0 * Math.atan(
        Math.sqrt((eccentricity + 1.0) / (eccentricity - 1.0)) * Math.tanh(hyperbolicAnomaly / 2.0)
    );
    
    const inclination = Math.PI * inclinationDeg / 180;
    const longAscNode = Math.PI * longAscNodeDeg / 180;
    const argPerihelion = Math.PI * argPerihelionDeg / 180;
    
    const cosO = Math.cos(longAscNode);
    const sinO = Math.sin(longAscNode);
    const cosI = Math.cos(inclination);
    const sinI = Math.sin(inclination);
    
    const argumentSum = argPerihelion + trueAnomaly;
    const cosArg = Math.cos(argumentSum);
    const sinArg = Math.sin(argumentSum);
    
    const x = radius * (cosO * cosArg - sinO * sinArg * cosI);
    const y = radius * (sinO * cosArg + cosO * sinArg * cosI);
    
    return [name, symbol, x, y];
}

/**
 * Compute planet positions using circular orbits (fallback)
 */
function computePlanetPositions(epoch) {
    const days = daysSinceEpoch(epoch);
    const positions = [];
    
    for (const [name, periodDays, radius, lon0, symbol] of PLANETS) {
        const meanMotion = 360.0 / periodDays;
        const angleDeg = (lon0 + meanMotion * days) % 360.0;
        const angle = Math.PI * angleDeg / 180;
        const x = Math.cos(angle) * radius;
        const y = Math.sin(angle) * radius;
        positions.push([name, symbol, x, y]);
    }
    
    const cometPosition = computeCometPosition(epoch);
    if (cometPosition) {
        positions.push(cometPosition);
    }
    
    return positions;
}

/**
 * Fetch planet position from NASA HORIZONS API
 */
async function fetchHorizonsPosition(planetName, epoch) {
    if (!(planetName in HORIZONS_PLANET_IDS)) return null;
    
    const planetId = HORIZONS_PLANET_IDS[planetName];
    let commandValue = String(planetId).trim();
    if (commandValue.includes(" ") && !/^-?\d+$/.test(commandValue)) {
        commandValue = `"${commandValue}"`;
    }
    
    const jd = datetimeToJulianDate(epoch);
    const jdStart = `JD${jd.toFixed(6)}`;
    const jdStop = `JD${(jd + 1/1440).toFixed(6)}`;
    
    const params = new URLSearchParams({
        format: "text",
        COMMAND: commandValue,
        OBJ_DATA: "NO",
        MAKE_EPHEM: "YES",
        EPHEM_TYPE: "VECTORS",
        CENTER: "500@10",
        START_TIME: jdStart,
        STOP_TIME: jdStop,
        STEP_SIZE: "1m",
        VEC_TABLE: "2",
        VEC_LABELS: "YES",
        CSV_FORMAT: "NO",
    });
    
    try {
        let url = `${HORIZONS_API_URL}?${params}`;
        
        // Use CORS proxy if enabled (required for browser requests)
        if (USE_CORS_PROXY) {
            url = `${CORS_PROXY_URL}${encodeURIComponent(url)}`;
        }
        
        const controller = new AbortController();
        const timeoutId = setTimeout(() => controller.abort(), HORIZONS_TIMEOUT);
        
        const response = await fetch(url, { signal: controller.signal });
        clearTimeout(timeoutId);
        
        if (!response.ok) {
            console.warn(`HORIZONS API fetch failed for ${planetName}: ${response.status} ${response.statusText}`);
            return null;
        }
        
        let resultText;
        if (USE_CORS_PROXY) {
            // CORS proxy returns JSON with contents field
            const data = await response.json();
            resultText = data.contents || "";
        } else {
            resultText = await response.text();
        }
        
        if (!resultText || resultText.includes("API Error") || resultText.includes("No matches found")) {
            console.warn(`HORIZONS API error for ${planetName}: ${resultText.substring(0, 200)}`);
            return null;
        }
        
        // Parse response
        const lines = resultText.split("\n");
        const KM_TO_AU = 1.0 / 149597870.7;
        let inData = false;
        
        for (let i = 0; i < lines.length; i++) {
            const line = lines[i];
            if (line.includes("$$SOE")) {
                inData = true;
                continue;
            }
            if (line.includes("$$EOE")) break;
            
            if (inData && line.trim()) {
                const lineUpper = line.toUpperCase();
                if ((lineUpper.includes("X=") || lineUpper.includes("X =")) && i + 1 < lines.length) {
                    const xMatch = lineUpper.match(/X\s*=\s*([+-]?\d+\.?\d*[Ee][+-]?\d+|[+-]?\d+\.?\d*)/);
                    const yMatch = lineUpper.match(/Y\s*=\s*([+-]?\d+\.?\d*[Ee][+-]?\d+|[+-]?\d+\.?\d*)/);
                    
                    if (xMatch && yMatch) {
                        const nextLineUpper = lines[i + 1].toUpperCase();
                        const vxMatch = nextLineUpper.match(/VX\s*=\s*([+-]?\d+\.?\d*[Ee][+-]?\d+|[+-]?\d+\.?\d*)/);
                        const vyMatch = nextLineUpper.match(/VY\s*=\s*([+-]?\d+\.?\d*[Ee][+-]?\d+|[+-]?\d+\.?\d*)/);
                        
                        if (vxMatch && vyMatch) {
                            try {
                                const xKm = parseFloat(xMatch[1]);
                                const yKm = parseFloat(yMatch[1]);
                                const vxKms = parseFloat(vxMatch[1]);
                                const vyKms = parseFloat(vyMatch[1]);
                                
                                const xAu = xKm * KM_TO_AU;
                                const yAu = yKm * KM_TO_AU;
                                const vxAuPerDay = (vxKms * 86400) * KM_TO_AU;
                                const vyAuPerDay = (vyKms * 86400) * KM_TO_AU;
                                
                                if (Math.abs(xAu) < 100 && Math.abs(yAu) < 100) {
                                    return [xAu, yAu, vxAuPerDay, vyAuPerDay];
                                }
                            } catch (e) {
                                continue;
                            }
                        }
                    }
                }
            }
        }
        
        return null;
    } catch (e) {
        return null;
    }
}

/**
 * Fetch all planet positions from HORIZONS API
 * Returns positions if at least the major planets are successfully fetched
 */
async function fetchPlanetPositionsOnline(epoch) {
    const positions = [];
    let successCount = 0;
    let failCount = 0;
    
    // Fetch major planets (require at least 4 to succeed)
    for (const [name, , , , symbol] of PLANETS) {
        const coords = await fetchHorizonsPosition(name, epoch);
        if (coords === null) {
            failCount++;
            console.warn(`Failed to fetch ${name} from HORIZONS API`);
            continue;
        }
        const [x, y, vx, vy] = coords;
        positions.push([name, symbol, x, y, vx, vy]);
        successCount++;
    }
    
    // Require at least 4 major planets to be successful
    if (successCount < 4) {
        console.warn(`Only ${successCount} planets fetched, need at least 4. Falling back to calculated.`);
        return null;
    }
    
    // Fetch additional objects (optional - continue even if they fail)
    for (const [name, symbol] of ADDITIONAL_INNER_OBJECTS) {
        const coords = await fetchHorizonsPosition(name, epoch);
        if (coords === null) continue;
        const [x, y, vx, vy] = coords;
        positions.push([name, symbol, x, y, vx, vy]);
    }
    
    for (const [name, symbol] of ADDITIONAL_OUTER_OBJECTS) {
        const coords = await fetchHorizonsPosition(name, epoch);
        if (coords === null) continue;
        const [x, y, vx, vy] = coords;
        positions.push([name, symbol, x, y, vx, vy]);
    }
    
    // Add comet (calculated, not from API)
    const cometPosition = computeCometPosition(epoch);
    if (cometPosition) {
        const [name, symbol, x, y] = cometPosition;
        positions.push([name, symbol, x, y, 0.0, 0.0]);
    }
    
    return positions;
}

/**
 * Interpolate positions based on velocity
 */
function interpolatePositions(basePositions, baseTime, currentTime) {
    const elapsedDays = (currentTime - baseTime) / (1000 * 60 * 60 * 24);
    const interpolated = [];
    
    for (const [name, symbol, xBase, yBase, vx, vy] of basePositions) {
        const xNew = xBase + vx * elapsedDays;
        const yNew = yBase + vy * elapsedDays;
        interpolated.push([name, symbol, xNew, yNew]);
    }
    
    return interpolated;
}

/**
 * Scale radius for logarithmic display
 */
function scaleRadius(radius, maxPixels, minRadius, maxRadius) {
    if (radius <= 0 || maxPixels <= 0) return 0.0;
    
    const minR = Math.max(minRadius, 1e-6);
    const maxR = Math.max(maxRadius, minR);
    const clampedRadius = Math.max(radius, minR);
    
    const logMin = Math.log(minR);
    const logMax = Math.log(maxR);
    const logRadius = Math.log(clampedRadius);
    
    let normalizedLog;
    if (Math.abs(logMin - logMax) < 1e-10) {
        normalizedLog = 1.0;
    } else {
        normalizedLog = (logRadius - logMin) / (logMax - logMin);
    }
    
    const normalized = MIN_RADIUS_FRACTION + (1.0 - MIN_RADIUS_FRACTION) * normalizedLog;
    return Math.max(MIN_RADIUS_FRACTION, Math.min(1.0, normalized)) * maxPixels;
}

/**
 * Compute radius bounds for a set of positions
 */
function computeRadiusBounds(positions) {
    if (positions.length === 0) {
        const radii = PLANETS.map(p => p[2]);
        return [Math.min(...radii), Math.max(...radii)];
    }
    
    const radii = positions.map(([, , x, y]) => Math.max(1e-6, Math.hypot(x, y)));
    let minRadius = Math.min(...radii);
    let maxRadius = Math.max(...radii);
    if (Math.abs(minRadius - maxRadius) < 1e-10) {
        maxRadius = maxRadius * 1.1;
    }
    return [minRadius, maxRadius];
}

/**
 * Render planet chart on canvas
 * All drawing is constrained to the canvas boundaries using clipping
 */
function renderPlanetChart(ctx, positions, minRadius, maxRadius) {
    if (!ctx) return;
    
    const width = ctx.canvas.width;
    const height = ctx.canvas.height;
    
    if (width < 6 || height < 4 || positions.length === 0) return;
    
    // Clear and fill background
    ctx.clearRect(0, 0, width, height);
    ctx.fillStyle = "#001100";
    ctx.fillRect(0, 0, width, height);
    
    // Set up clipping region to constrain all drawing to canvas boundaries
    ctx.save();
    ctx.beginPath();
    ctx.rect(0, 0, width, height);
    ctx.clip();
    
    const centerX = width / 2;
    const centerY = height / 2;
    
    // Calculate maximum radius that fits within canvas bounds
    // Account for both horizontal and vertical constraints
    const maxHorizontalRadius = (width / 2) - 2; // Leave 2px margin
    const maxVerticalRadius = ((height / 2) / CHAR_ASPECT_RATIO) - 2; // Account for aspect ratio
    const maxRadiusPx = Math.max(1, Math.min(maxHorizontalRadius, maxVerticalRadius));
    
    // Draw orbit paths (will be clipped automatically)
    const drawnPaths = new Set();
    for (const [name, symbol, x, y] of positions) {
        const radius = Math.hypot(x, y);
        const radiusPx = scaleRadius(radius, maxRadiusPx, minRadius, maxRadius);
        const radiusInt = Math.round(radiusPx);
        
        // Only draw if radius fits within bounds
        if (radiusInt > 1 && radiusPx <= maxRadiusPx && !drawnPaths.has(radiusInt)) {
            ctx.strokeStyle = "#003300";
            ctx.lineWidth = 1;
            ctx.beginPath();
            ctx.arc(centerX, centerY, radiusPx, 0, 2 * Math.PI);
            ctx.stroke();
            drawnPaths.add(radiusInt);
        }
    }
    
    // Draw sun (ensure it's within bounds)
    const sunRadius = Math.min(3, maxRadiusPx - 1);
    if (sunRadius > 0 && centerX - sunRadius >= 0 && centerX + sunRadius < width &&
        centerY - sunRadius >= 0 && centerY + sunRadius < height) {
        ctx.fillStyle = "#ffff00";
        ctx.beginPath();
        ctx.arc(centerX, centerY, sunRadius, 0, 2 * Math.PI);
        ctx.fill();
        ctx.strokeStyle = "#ffff00";
        ctx.lineWidth = 1;
        const crossSize = Math.min(5, maxRadiusPx - 1);
        if (crossSize > 0) {
            ctx.beginPath();
            ctx.moveTo(Math.max(0, centerX - crossSize), centerY);
            ctx.lineTo(Math.min(width, centerX + crossSize), centerY);
            ctx.moveTo(centerX, Math.max(0, centerY - crossSize));
            ctx.lineTo(centerX, Math.min(height, centerY + crossSize));
            ctx.stroke();
        }
    }
    
    // Draw planets (all constrained by clipping)
    const occupied = new Set();
    for (const [name, symbol, x, y] of positions) {
        const radius = Math.hypot(x, y);
        const radiusPx = scaleRadius(radius, maxRadiusPx, minRadius, maxRadius);
        const angle = Math.atan2(y, x);
        const drawX = centerX + Math.cos(angle) * radiusPx;
        const drawY = centerY - Math.sin(angle) * radiusPx * CHAR_ASPECT_RATIO;
        
        // Clamp coordinates to canvas bounds
        const clampedX = Math.max(0, Math.min(width - 1, drawX));
        const clampedY = Math.max(0, Math.min(height - 1, drawY));
        
        // Find available position within bounds
        let placed = false;
        for (let dy = -2; dy <= 2 && !placed; dy++) {
            for (let dx = -2; dx <= 2 && !placed; dx++) {
                const px = Math.round(clampedX + dx);
                const py = Math.round(clampedY + dy);
                const key = `${py},${px}`;
                
                // Ensure position is within canvas bounds
                if (px >= 0 && px < width && py >= 0 && py < height && !occupied.has(key)) {
                    ctx.fillStyle = PLANET_COLORS[name] || "#00ff00";
                    ctx.font = "12px monospace";
                    ctx.textAlign = "center";
                    ctx.textBaseline = "middle";
                    ctx.fillText(symbol, px, py);
                    occupied.add(key);
                    placed = true;
                }
            }
        }
        
        // Fallback: draw at clamped position if not placed
        if (!placed) {
            const px = Math.round(clampedX);
            const py = Math.round(clampedY);
            if (px >= 0 && px < width && py >= 0 && py < height) {
                ctx.fillStyle = PLANET_COLORS[name] || "#00ff00";
                ctx.font = "12px monospace";
                ctx.textAlign = "center";
                ctx.textBaseline = "middle";
                ctx.fillText(symbol, px, py);
            }
        }
    }
    
    // Restore clipping (remove clipping region)
    ctx.restore();
}

/**
 * Fetch weather data from Open-Meteo API
 */
async function fetchWeatherData(lat, lon) {
    const params = new URLSearchParams({
        latitude: lat.toFixed(4),
        longitude: lon.toFixed(4),
        current_weather: "true",
        hourly: "relativehumidity_2m,apparent_temperature,dewpoint_2m,precipitation,rain,snowfall,cloudcover,visibility,surface_pressure,pressure_msl",
        daily: "sunrise,sunset,uv_index_max,precipitation_probability_max",
        temperature_unit: "celsius",
        windspeed_unit: "kmh",
        timezone: "UTC",
        forecast_days: "1",
    });
    
    try {
        const url = `${WEATHER_API_URL}?${params}`;
        const response = await fetch(url);
        if (!response.ok) return null;
        
        const data = await response.json();
        const current = data.current_weather || {};
        if (!current.time) return null;
        
        const timestamp = new Date(current.time);
        const hourly = data.hourly || {};
        const hourlyTimes = hourly.time || [];
        const hourlyIndex = hourlyTimes.indexOf(current.time);
        
        const extractHourly = (key) => {
            if (hourlyIndex === -1) return null;
            const series = hourly[key];
            if (!Array.isArray(series) || hourlyIndex >= series.length) return null;
            return typeof series[hourlyIndex] === 'number' ? series[hourlyIndex] : null;
        };
        
        const daily = data.daily || {};
        const dailyTimes = daily.time || [];
        const dailyIndex = dailyTimes.length > 0 ? dailyTimes.indexOf(timestamp.toISOString().split('T')[0]) : 0;
        
        const extractDaily = (key) => {
            if (dailyIndex === -1) return null;
            const series = daily[key];
            if (!Array.isArray(series) || dailyIndex >= series.length) return null;
            const value = series[dailyIndex];
            if (typeof value === 'string') {
                try {
                    return new Date(value);
                } catch {
                    return null;
                }
            }
            return typeof value === 'number' ? value : null;
        };
        
        const tempC = current.temperature;
        const tempF = tempC * 9 / 5 + 32;
        const apparentC = extractHourly('apparent_temperature') ?? tempC;
        const apparentF = apparentC * 9 / 5 + 32;
        const windKmh = current.windspeed;
        const windMph = windKmh * 0.621371;
        const windDirDeg = current.winddirection;
        const windCardinal = windDirectionFromDegrees(windDirDeg);
        const humidity = extractHourly('relativehumidity_2m');
        const dewpointC = extractHourly('dewpoint_2m');
        const dewpointF = dewpointC !== null ? dewpointC * 9 / 5 + 32 : null;
        const cloudCover = extractHourly('cloudcover');
        const visibilityM = extractHourly('visibility');
        const visibilityKm = visibilityM !== null ? visibilityM / 1000.0 : null;
        const visibilityMiles = visibilityKm !== null ? visibilityKm * 0.621371 : null;
        const surfacePressureHpa = extractHourly('surface_pressure');
        const surfacePressureInhg = surfacePressureHpa !== null ? surfacePressureHpa * 0.0295299830714 : null;
        const pressureMslHpa = extractHourly('pressure_msl');
        const pressureMslInhg = pressureMslHpa !== null ? pressureMslHpa * 0.0295299830714 : null;
        const precipMm = extractHourly('precipitation');
        const precipIn = precipMm !== null ? precipMm * 0.0393701 : null;
        const rainMm = extractHourly('rain');
        const rainIn = rainMm !== null ? rainMm * 0.0393701 : null;
        const snowMm = extractHourly('snowfall');
        const snowIn = snowMm !== null ? snowMm * 0.0393701 : null;
        const uvIndexMax = extractDaily('uv_index_max');
        const precipProbabilityMax = extractDaily('precipitation_probability_max');
        const sunrise = extractDaily('sunrise');
        const sunset = extractDaily('sunset');
        const isDay = current.is_day === 1;
        const weatherCode = current.weathercode;
        const condition = WEATHER_CODE_MAP[weatherCode] || "Unknown";
        
        return {
            timestamp,
            temperature_c: tempC,
            temperature_f: tempF,
            apparent_c: apparentC,
            apparent_f: apparentF,
            windspeed_kmh: windKmh,
            windspeed_mph: windMph,
            wind_direction_deg: windDirDeg,
            wind_direction_cardinal: windCardinal,
            humidity,
            dewpoint_c: dewpointC,
            dewpoint_f: dewpointF,
            cloud_cover: cloudCover,
            visibility_km: visibilityKm,
            visibility_miles: visibilityMiles,
            surface_pressure_hpa: surfacePressureHpa,
            surface_pressure_inhg: surfacePressureInhg,
            pressure_msl_hpa: pressureMslHpa,
            pressure_msl_inhg: pressureMslInhg,
            precipitation_mm: precipMm,
            precipitation_in: precipIn,
            rain_mm: rainMm,
            rain_in: rainIn,
            snow_mm: snowMm,
            snow_in: snowIn,
            uv_index_max: uvIndexMax,
            precipitation_probability_max: precipProbabilityMax,
            sunrise,
            sunset,
            is_day: isDay,
            condition,
            weather_code: weatherCode,
        };
    } catch (e) {
        return null;
    }
}

/**
 * Convert wind direction degrees to cardinal direction
 */
function windDirectionFromDegrees(degrees) {
    if (degrees === null || degrees === undefined) return "--";
    const bearings = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"];
    const index = Math.floor(((degrees % 360) + 22.5) / 45) % bearings.length;
    return bearings[index];
}

/**
 * Load weather coordinates with priority: GET params > config.js defaults > localStorage
 * 
 * Priority order:
 * 1. URL GET parameters (?lat=X&lon=Y) - highest priority, allows runtime override
 * 2. config.js defaults (from .env file) - default values
 * 3. localStorage - user preference persistence
 */
function loadWeatherCoordinates() {
    const urlParams = new URLSearchParams(window.location.search);
    
    // Priority 1: GET parameters (highest priority - runtime override)
    let lat = urlParams.get('lat');
    let lon = urlParams.get('lon');
    
    // Priority 2: config.js defaults (from .env file)
    if (!lat && typeof DEFAULT_WEATHER_LAT !== 'undefined' && DEFAULT_WEATHER_LAT !== null) {
        lat = String(DEFAULT_WEATHER_LAT);
    }
    if (!lon && typeof DEFAULT_WEATHER_LON !== 'undefined' && DEFAULT_WEATHER_LON !== null) {
        lon = String(DEFAULT_WEATHER_LON);
    }
    
    // Priority 3: localStorage (user preference)
    if (!lat) {
        lat = localStorage.getItem(WEATHER_LAT_ENV);
    }
    if (!lon) {
        lon = localStorage.getItem(WEATHER_LON_ENV);
    }
    
    if (lat && lon) {
        const latNum = parseFloat(lat);
        const lonNum = parseFloat(lon);
        if (!isNaN(latNum) && !isNaN(lonNum)) {
            return [latNum, lonNum];
        }
    }
    return null;
}

/**
 * Fetch latest headlines from RSS feeds
 */
async function fetchLatestHeadlines() {
    const aggregated = [];
    
    for (const [sourceName, feedUrl] of NEWS_FEEDS) {
        try {
            // Use CORS proxy for RSS feeds (in production, use your own proxy)
            const proxyUrl = `https://api.allorigins.win/get?url=${encodeURIComponent(feedUrl)}`;
            const response = await fetch(proxyUrl);
            if (!response.ok) continue;
            
            const data = await response.json();
            const text = data.contents;
            
            const parser = new DOMParser();
            const doc = parser.parseFromString(text, 'text/xml');
            const items = doc.querySelectorAll('item');
            
            for (const item of items) {
                const title = item.querySelector('title')?.textContent?.trim();
                if (!title) continue;
                
                const pubDate = item.querySelector('pubDate')?.textContent;
                let published = null;
                let timestamp = "----";
                
                if (pubDate) {
                    try {
                        published = new Date(pubDate);
                        timestamp = published.toISOString().substr(11, 5) + 'Z';
                    } catch (e) {
                        // Keep default timestamp
                    }
                }
                
                aggregated.push({
                    published: published || new Date(0),
                    text: `[${sourceName} ${timestamp}] ${title}`
                });
            }
        } catch (e) {
            continue;
        }
    }
    
    if (aggregated.length === 0) return [];
    
    aggregated.sort((a, b) => a.published - b.published);
    if (aggregated.length > NEWS_MAX_ITEMS) {
        aggregated.splice(0, aggregated.length - NEWS_MAX_ITEMS);
    }
    
    // Shuffle
    for (let i = aggregated.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [aggregated[i], aggregated[j]] = [aggregated[j], aggregated[i]];
    }
    
    return aggregated.slice(0, NEWS_MAX_ITEMS).map(item => item.text);
}

/**
 * Update power bar buffer
 */
function updatePowerBarBuffer(buffer, width) {
    if (width <= 0) return "";
    
    if (buffer.length < width) {
        buffer = (buffer + " ".repeat(width)).slice(-width);
    } else if (buffer.length > width) {
        buffer = buffer.slice(-width);
    }
    
    const chars = buffer.split('');
    chars.shift();
    chars.push(POWER_BAR_TEXTURES[Math.floor(Math.random() * POWER_BAR_TEXTURES.length)]);
    
    for (let i = 0; i < Math.max(1, Math.floor(width / 12)); i++) {
        const idx = Math.floor(Math.random() * width);
        chars[idx] = POWER_BAR_TEXTURES[Math.floor(Math.random() * POWER_BAR_TEXTURES.length)];
    }
    
    return chars.join('');
}

/**
 * Render weather panel
 */
function renderWeatherPanel(weatherInfo, weatherError, animationPhase, powerBarState) {
    const content = document.getElementById('weather-content');
    const powerSection = document.getElementById('power-bar-section');
    const powerStatus = document.getElementById('power-status');
    const powerBar = document.getElementById('power-bar');
    const powerOutput = document.getElementById('power-output');
    
    if (!content) return powerBarState;
    
    content.innerHTML = '';
    
    if (weatherInfo) {
        const lines = [];
        lines.push(`Conditions: ${weatherInfo.condition || '--'}`);
        
        if (weatherInfo.dewpoint_c !== null && weatherInfo.dewpoint_f !== null) {
            lines.push(`Dew point: ${weatherInfo.dewpoint_c.toFixed(1)}°C (${weatherInfo.dewpoint_f.toFixed(1)}°F)`);
        }
        
        if (weatherInfo.cloud_cover !== null) {
            lines.push(`Clouds: ${Math.round(weatherInfo.cloud_cover)}% cover`);
        }
        
        if (weatherInfo.visibility_km !== null) {
            if (weatherInfo.visibility_miles !== null) {
                lines.push(`Visibility: ${weatherInfo.visibility_km.toFixed(1)} km (${weatherInfo.visibility_miles.toFixed(1)} mi)`);
            } else {
                lines.push(`Visibility: ${weatherInfo.visibility_km.toFixed(1)} km`);
            }
        }
        
        if (weatherInfo.surface_pressure_hpa !== null) {
            let pressureLine = `Pressure: ${Math.round(weatherInfo.surface_pressure_hpa)} hPa`;
            if (weatherInfo.surface_pressure_inhg !== null) {
                pressureLine += ` (${weatherInfo.surface_pressure_inhg.toFixed(2)} inHg)`;
            }
            lines.push(pressureLine);
        }
        
        if (weatherInfo.pressure_msl_hpa !== null) {
            let mslLine = `MSL Pressure: ${Math.round(weatherInfo.pressure_msl_hpa)} hPa`;
            if (weatherInfo.pressure_msl_inhg !== null) {
                mslLine += ` (${weatherInfo.pressure_msl_inhg.toFixed(2)} inHg)`;
            }
            lines.push(mslLine);
        }
        
        if (weatherInfo.precipitation_mm !== null) {
            let precipLine = `Precip (1h): ${weatherInfo.precipitation_mm.toFixed(1)} mm`;
            if (weatherInfo.precipitation_in !== null) {
                precipLine += ` (${weatherInfo.precipitation_in.toFixed(2)} in)`;
            }
            lines.push(precipLine);
        }
        
        if (weatherInfo.rain_mm !== null && weatherInfo.rain_mm > 0) {
            let rainLine = `Rain (1h): ${weatherInfo.rain_mm.toFixed(1)} mm`;
            if (weatherInfo.rain_in !== null) {
                rainLine += ` (${weatherInfo.rain_in.toFixed(2)} in)`;
            }
            lines.push(rainLine);
        }
        
        if (weatherInfo.snow_mm !== null && weatherInfo.snow_mm > 0) {
            let snowLine = `Snow (1h): ${weatherInfo.snow_mm.toFixed(1)} mm`;
            if (weatherInfo.snow_in !== null) {
                snowLine += ` (${weatherInfo.snow_in.toFixed(2)} in)`;
            }
            lines.push(snowLine);
        }
        
        if (weatherInfo.sunrise instanceof Date) {
            lines.push(`Sunrise: ${weatherInfo.sunrise.toISOString().substr(11, 5)} UTC`);
        }
        
        if (weatherInfo.sunset instanceof Date) {
            lines.push(`Sunset: ${weatherInfo.sunset.toISOString().substr(11, 5)} UTC`);
        }
        
        if (typeof weatherInfo.is_day === 'boolean') {
            lines.push(`Daylight: ${weatherInfo.is_day ? 'Yes' : 'No'}`);
        }
        
        // Gauge lines
        const gaugeLines = [];
        
        if (weatherInfo.temperature_c !== null) {
            const tempDisplay = weatherInfo.temperature_f !== null
                ? `${weatherInfo.temperature_c.toFixed(1)}°C / ${weatherInfo.temperature_f.toFixed(1)}°F`
                : `${weatherInfo.temperature_c.toFixed(1)}°C`;
            gaugeLines.push(`Temp: ${tempDisplay}`);
        }
        
        if (weatherInfo.apparent_c !== null) {
            const feelsDisplay = weatherInfo.apparent_f !== null
                ? `${weatherInfo.apparent_c.toFixed(1)}°C / ${weatherInfo.apparent_f.toFixed(1)}°F`
                : `${weatherInfo.apparent_c.toFixed(1)}°C`;
            gaugeLines.push(`Feels: ${feelsDisplay}`);
        }
        
        if (weatherInfo.windspeed_mph !== null || weatherInfo.windspeed_kmh !== null) {
            const windParts = [];
            if (weatherInfo.windspeed_mph !== null) windParts.push(`${weatherInfo.windspeed_mph.toFixed(1)} mph`);
            if (weatherInfo.windspeed_kmh !== null) windParts.push(`${weatherInfo.windspeed_kmh.toFixed(1)} km/h`);
            let windDisplay = windParts.join(' / ');
            if (weatherInfo.wind_direction_cardinal && weatherInfo.wind_direction_cardinal !== '--') {
                windDisplay += ` ${weatherInfo.wind_direction_cardinal}`;
            }
            gaugeLines.push(`Wind: ${windDisplay}`);
        }
        
        if (weatherInfo.humidity !== null) {
            gaugeLines.push(`Humidity: ${Math.round(weatherInfo.humidity)}%`);
        }
        
        if (weatherInfo.precipitation_probability_max !== null) {
            gaugeLines.push(`Precip prob: ${Math.round(weatherInfo.precipitation_probability_max)}%`);
        }
        
        if (weatherInfo.uv_index_max !== null) {
            gaugeLines.push(`UV: ${weatherInfo.uv_index_max.toFixed(1)}`);
        }
        
        if (gaugeLines.length > 0 && lines.length > 0) {
            lines.unshift('');
        }
        
        lines.unshift(...gaugeLines);
        
        for (const line of lines) {
            const div = document.createElement('div');
            div.className = 'weather-line';
            div.textContent = line;
            content.appendChild(div);
        }
    } else {
        const message = weatherError || `Set ${WEATHER_LAT_ENV}/${WEATHER_LON_ENV} URL params or localStorage to enable weather.`;
        const div = document.createElement('div');
        div.className = 'weather-message';
        div.textContent = message;
        content.appendChild(div);
    }
    
    // Update power bar (only update spinner, not the bar itself - bar updates separately)
    if (powerSection && powerStatus) {
        const spinner = WEATHER_SPINNER_CHARS[animationPhase % WEATHER_SPINNER_CHARS.length];
        powerStatus.textContent = `Power Feed ${spinner}`;
    }
    
    return powerBarState;
}

/**
 * Render telemetry panel
 */
function renderTelemetryPanel(elementId, positions, currentTime, highlightDuration) {
    const content = document.getElementById(elementId);
    if (!content) return;
    
    content.innerHTML = '';
    
    if (positions.length === 0) {
        const div = document.createElement('div');
        div.className = 'telemetry-message';
        div.textContent = 'No telemetry data available';
        content.appendChild(div);
        return;
    }
    
    // Sort by distance
    const positionsWithDistance = positions.map(([name, symbol, x, y]) => {
        const distance = Math.hypot(x, y);
        return [name, symbol, x, y, distance];
    });
    positionsWithDistance.sort((a, b) => a[4] - b[4]);
    
    // Update change tracking
    if (currentTime !== null) {
        for (const [name, , x, y, distance] of positionsWithDistance) {
            const key = name;
            if (key in previousTelemetry) {
                const prev = previousTelemetry[key];
                const xChanged = Math.abs(x - prev.x) > 1e-8;
                const yChanged = Math.abs(y - prev.y) > 1e-8;
                const distChanged = Math.abs(distance - prev.dist) > 1e-8;
                
                if (xChanged) telemetryChangeTimes[`${key}_x`] = currentTime;
                if (yChanged) telemetryChangeTimes[`${key}_y`] = currentTime;
                if (distChanged) telemetryChangeTimes[`${key}_dist`] = currentTime;
                
                previousTelemetry[key] = { x, y, dist: distance, timestamp: currentTime };
            } else {
                previousTelemetry[key] = { x, y, dist: distance, timestamp: currentTime };
                telemetryChangeTimes[`${key}_x`] = currentTime;
                telemetryChangeTimes[`${key}_y`] = currentTime;
                telemetryChangeTimes[`${key}_dist`] = currentTime;
            }
        }
    }
    
    // Render lines
    for (const [name, symbol, x, y, distance] of positionsWithDistance) {
        const xStr = x >= 0 ? `+${x.toFixed(6)}` : x.toFixed(6);
        const yStr = y >= 0 ? `+${y.toFixed(6)}` : y.toFixed(6);
        const distStr = distance.toFixed(6);
        
        const key = name;
        const xKey = `${key}_x`;
        const yKey = `${key}_y`;
        const distKey = `${key}_dist`;
        
        let xHighlight = false;
        let yHighlight = false;
        let distHighlight = false;
        
        if (currentTime !== null) {
            if (xKey in telemetryChangeTimes) {
                const timeSinceChange = currentTime - telemetryChangeTimes[xKey];
                xHighlight = timeSinceChange < highlightDuration && timeSinceChange >= 0;
            }
            if (yKey in telemetryChangeTimes) {
                const timeSinceChange = currentTime - telemetryChangeTimes[yKey];
                yHighlight = timeSinceChange < highlightDuration && timeSinceChange >= 0;
            }
            if (distKey in telemetryChangeTimes) {
                const timeSinceChange = currentTime - telemetryChangeTimes[distKey];
                distHighlight = timeSinceChange < highlightDuration && timeSinceChange >= 0;
            }
        }
        
        const prefix = `${symbol} ${name.padEnd(8).substring(0, 8)} `;
        const xPart = `X:${xStr.padStart(12)} `;
        const yPart = `Y:${yStr.padStart(12)} `;
        const distPart = `D:${distStr.padStart(10)}AU`;
        
        const line = document.createElement('div');
        line.className = 'telemetry-line';
        
        const prefixSpan = document.createElement('span');
        prefixSpan.textContent = prefix;
        line.appendChild(prefixSpan);
        
        const xSpan = document.createElement('span');
        xSpan.textContent = xPart;
        if (xHighlight) xSpan.className = 'highlight';
        line.appendChild(xSpan);
        
        const ySpan = document.createElement('span');
        ySpan.textContent = yPart;
        if (yHighlight) ySpan.className = 'highlight';
        line.appendChild(ySpan);
        
        const distSpan = document.createElement('span');
        distSpan.textContent = distPart;
        if (distHighlight) distSpan.className = 'highlight';
        line.appendChild(distSpan);
        
        content.appendChild(line);
    }
}

/**
 * Render news panel
 */
function renderNewsPanel(newsHistory, newsPartialLine) {
    const content = document.getElementById('news-content');
    if (!content) return;
    
    content.innerHTML = '';
    
    const allLines = [...newsHistory];
    if (newsPartialLine) {
        allLines.push(newsPartialLine + '<span class="news-cursor">_</span>');
    }
    
    if (allLines.length === 0) {
        const div = document.createElement('div');
        div.className = 'news-message';
        div.textContent = 'Awaiting latest headlines...';
        content.appendChild(div);
        return;
    }
    
    for (const line of allLines) {
        const div = document.createElement('div');
        div.className = 'news-item';
        if (line.includes('<span')) {
            div.innerHTML = line;
            div.classList.add('news-current');
        } else {
            div.textContent = line;
        }
        content.appendChild(div);
    }
    
    // Scroll to bottom
    content.scrollTop = content.scrollHeight;
}

/**
 * Render inner chart
 */
function renderInnerChart(positions, lastRefresh, secondsSinceRefresh, dataSource) {
    const header = document.getElementById('inner-chart-header');
    if (header) {
        const refreshText = lastRefresh
            ? `Last refresh ${lastRefresh.toISOString().substr(11, 8)} UTC`
            : 'Waiting for first update...';
        const nextUpdate = Math.max(0, SYNC_UPDATE_INTERVAL_SECONDS - secondsSinceRefresh);
        header.textContent = `Solar System (top-down) - ${refreshText} - Next sync update in ${Math.floor(nextUpdate)}s - [${dataSource}]`;
    }
    
    const innerPositions = positions.filter(([name]) => INNER_PLANET_NAMES.has(name));
    if (innerPositions.length === 0) return;
    
    const [minRadius, maxRadius] = computeRadiusBounds(innerPositions);
    renderPlanetChart(innerChartCtx, innerPositions, minRadius * 0.95, maxRadius * 1.25);
}

/**
 * Render outer chart
 */
function renderOuterChart(positions) {
    const outerPositions = positions.filter(([name]) => OUTER_PLANET_NAMES.has(name));
    if (outerPositions.length === 0) return;
    
    const [minRadius, maxRadius] = computeRadiusBounds(outerPositions);
    renderPlanetChart(outerChartCtx, outerPositions, minRadius * 0.9, maxRadius * 1.1);
}

/**
 * Main render function
 */
function renderDashboard() {
    const now = new Date();
    const nowTs = now.getTime();
    const elapsed = (nowTs - lastSyncUpdateTs) / 1000;
    
    // Render charts
    renderInnerChart(positions, lastRefreshDt, elapsed, dataSource);
    renderOuterChart(positions);
    
    // Render weather
    powerBarState = renderWeatherPanel(weatherData, weatherError, weatherAnimationPhase, powerBarState);
    
    // Render news
    const newsPartial = newsCurrent ? newsCurrent.substring(0, newsProgress) + '_' : '';
    renderNewsPanel(newsHistory, newsPartial);
    
    // Render telemetry
    const innerPositions = positions.filter(([name]) => INNER_PLANET_NAMES.has(name));
    const outerPositions = positions.filter(([name]) => OUTER_PLANET_NAMES.has(name));
    const extras = positions.filter(([name]) => !INNER_PLANET_NAMES.has(name) && !OUTER_PLANET_NAMES.has(name));
    const allOuterPositions = [...outerPositions, ...extras];
    
    renderTelemetryPanel('inner-telemetry-content', innerPositions, nowTs / 1000, TELEMETRY_HIGHLIGHT_DURATION);
    renderTelemetryPanel('outer-telemetry-content', allOuterPositions, nowTs / 1000, TELEMETRY_HIGHLIGHT_DURATION);
}

/**
 * Main update loop
 */
async function updateData() {
    const now = new Date();
    const nowTs = now.getTime();
    const timeSinceSync = (nowTs - lastSyncUpdateTs) / 1000;
    
    // Synchronized update every 5 minutes
    if (positions.length === 0 || timeSinceSync >= SYNC_UPDATE_INTERVAL_SECONDS) {
        lastRefreshDt = now;
        lastSyncUpdateTs = nowTs;
        
        // Update planet positions
        try {
            const onlinePositions = await fetchPlanetPositionsOnline(lastRefreshDt);
            if (onlinePositions !== null && onlinePositions.length > 0) {
                basePositionsWithVelocities = onlinePositions;
                basePositionsTimestamp = lastRefreshDt;
                positions = interpolatePositions(onlinePositions, lastRefreshDt, now);
                dataSource = "NASA HORIZONS";
                console.log(`✓ Successfully loaded ${onlinePositions.length} positions from NASA HORIZONS API`);
            } else {
                positions = computePlanetPositions(lastRefreshDt);
                basePositionsWithVelocities = [];
                basePositionsTimestamp = null;
                dataSource = "calculated";
                console.warn("⚠ Falling back to calculated positions (HORIZONS API unavailable)");
            }
        } catch (e) {
            console.error("Error fetching HORIZONS data:", e);
            positions = computePlanetPositions(lastRefreshDt);
            basePositionsWithVelocities = [];
            basePositionsTimestamp = null;
            dataSource = "calculated";
        }
        
        // Update weather
        const weatherCoords = loadWeatherCoordinates();
        if (weatherCoords) {
            const [lat, lon] = weatherCoords;
            const fetched = await fetchWeatherData(lat, lon);
            if (fetched) {
                weatherData = fetched;
                weatherError = null;
            } else {
                weatherError = "Weather lookup failed (network or API issue).";
            }
        } else {
            weatherData = null;
            weatherError = `Configure ${WEATHER_LAT_ENV}/${WEATHER_LON_ENV} URL params or localStorage to show local weather.`;
        }
        
        // Update news
        const fetched = await fetchLatestHeadlines();
        if (fetched.length > 0) {
            newsItems = fetched;
            newsCurrent = "";
            newsProgress = 0;
            newsLastUpdateDt = now;
        }
    }
    
    // Interpolate positions continuously
    if (basePositionsWithVelocities.length > 0 && basePositionsTimestamp) {
        positions = interpolatePositions(basePositionsWithVelocities, basePositionsTimestamp, now);
    } else if (positions.length === 0 && lastRefreshDt) {
        positions = computePlanetPositions(lastRefreshDt);
    }
    
    // News ticker animation
    if (newsItems.length > 0 && !newsCurrent) {
        const index = Math.floor((nowTs / 1000) % newsItems.length);
        newsCurrent = newsItems[index];
        newsProgress = 0;
        newsCharTs = nowTs;
    }
    
    if (newsCurrent && nowTs - newsCharTs >= NEWS_CHAR_INTERVAL * 1000) {
        const elapsedChars = Math.max(1, Math.floor((nowTs - newsCharTs) / (NEWS_CHAR_INTERVAL * 1000)));
        newsProgress = Math.min(newsCurrent.length, newsProgress + elapsedChars);
        newsCharTs = nowTs;
        if (newsProgress >= newsCurrent.length) {
            newsHistory.push(newsCurrent);
            if (newsHistory.length > NEWS_HISTORY_LIMIT) {
                newsHistory = newsHistory.slice(-NEWS_HISTORY_LIMIT);
            }
            newsCurrent = "";
            newsProgress = 0;
        }
    }
    
    weatherAnimationPhase = (weatherAnimationPhase + 1) % (WEATHER_SPINNER_CHARS.length * 64);
}

/**
 * Update power bar separately at a slower rate
 */
function updatePowerBar() {
    const powerBar = document.getElementById('power-bar');
    const powerOutput = document.getElementById('power-output');
    const weatherContent = document.getElementById('weather-content');
    
    if (powerBar && powerOutput && weatherContent) {
        const barWidth = Math.max(8, weatherContent.offsetWidth - 20);
        powerBarState = updatePowerBarBuffer(powerBarState, barWidth);
        powerBar.textContent = `[${powerBarState}]`;
        
        const activeCells = powerBarState.split('').filter(ch => ch !== ' ' && ch !== '.').length;
        const utilization = powerBarState.length > 0 ? Math.floor((activeCells / powerBarState.length) * 100) : 0;
        powerOutput.textContent = `Output: ${utilization.toString().padStart(3, ' ')}%`;
    }
}

/**
 * Main animation loop
 */
let lastFrameTime = 0;
function mainLoop(currentTime) {
    const now = new Date();
    const nowTs = now.getTime();
    
    // Update data (async, but don't wait)
    if (nowTs - lastFrameTime >= 100) { // Update data every 100ms
        updateData();
        lastFrameTime = nowTs;
    }
    
    // Update power bar at slower rate
    if (nowTs - lastPowerBarUpdate >= POWER_BAR_UPDATE_INTERVAL) {
        updatePowerBar();
        lastPowerBarUpdate = nowTs;
    }
    
    // Render every frame
    renderDashboard();
    
    requestAnimationFrame(mainLoop);
}

// Initialize on page load
document.addEventListener('DOMContentLoaded', initDashboard);

// Handle keyboard shortcuts
document.addEventListener('keydown', (e) => {
    if (e.key === 'q' || e.key === 'Q') {
        // Could add quit functionality if needed
    }
    if (e.key === 'r' || e.key === 'R') {
        // Force refresh
        lastSyncUpdateTs = 0;
    }
});

