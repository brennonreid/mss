// asteroids.js
// Uses your customtrig.js Taylor-based unit vectors (angles in radians, no fixed-point).

// ---------- Game Constants ----------
const CANVAS_WIDTH = 800;
const CANVAS_HEIGHT = 600;
const SHIP_SIZE = 30; // diameter for drawing

// Rotation speed (~200 deg/sec) in radians/sec using TAU from customtrig
const ROTATION_RAD_PER_SEC = (200 / 360) * TAU;

const BULLET_SPEED_PIXELS_PER_SEC = 300;
const BULLET_LIFESPAN_SECONDS = 3;

const ASTEROID_SPEED_PIXELS_PER_SEC = 80;
const ASTEROID_VERTICES = 10;
const ASTEROID_JAGGEDNESS = 0.4;
const ASTEROID_MAX_SIZE = 60;
const ASTEROID_MIN_SIZE = 20;
const INITIAL_ASTEROIDS = 4;

// Max asteroid rotation speed (~60 deg/sec) in radians/sec
const ASTEROID_MAX_ROTATION_RAD_PER_SEC = (60 / 360) * TAU;

// Fixed timestep (60 Hz)
const PHYSICS_STEP = 1000 / 60;

// ---------- Game State ----------
let canvas, ctx;

const ship = {
  x: CANVAS_WIDTH / 2,
  y: CANVAS_HEIGHT / 2,
  angle: 0,                 // radians, 0 = facing +X
  rotationSpeed: 0,         // radians/sec (+/- ROTATION_RAD_PER_SEC)
  velocity: { x: 0, y: 0 },
  thrust_acceleration: 100, // px/s^2
  friction_coefficient: 0.1,
  canShoot: true,
  shootCooldown_seconds: 0.25,
  currentShootCooldown_seconds: 0,

  // interpolation
  prev_x: CANVAS_WIDTH / 2,
  prev_y: CANVAS_HEIGHT / 2,
  prev_angle: 0,
};

const bullets = [];
const asteroids = [];

// ---------- Timing ----------
let lastTime = 0;
let accumulatedTime = 0;

// ---------- Input ----------
const keys = { ArrowLeft: false, ArrowRight: false, ArrowUp: false, Space: false };

document.addEventListener('keydown', (e) => {
  if (e.code in keys) keys[e.code] = true;
});

document.addEventListener('keyup', (e) => {
  if (e.code in keys) keys[e.code] = false;
});

// ---------- Utilities ----------
function normalizeAngle(a) {
  a %= TAU;
  return a < 0 ? a + TAU : a;
}

function shortestAngleDelta(a, b) {
  let d = (b - a) % TAU;
  if (d >  Math.PI) d -= TAU;
  if (d < -Math.PI) d += TAU;
  return d;
}

function dirFromAngle(angle, precise = false) {
  // getUnitVectorFromAngle2 uses your Taylor + anchor system
  return getUnitVectorFromAngle2(angle, precise); // {x, y}
}

// ---------- Init ----------
window.onload = () => {
  canvas = document.getElementById('gameCanvas');
  ctx = canvas.getContext('2d');
  canvas.width = CANVAS_WIDTH;
  canvas.height = CANVAS_HEIGHT;

  for (let i = 0; i < INITIAL_ASTEROIDS; i++) spawnAsteroid(ASTEROID_MAX_SIZE);

  requestAnimationFrame(gameLoop);
};

// ---------- Game Loop ----------
function gameLoop(currentTime) {
  const deltaTime = currentTime - lastTime;
  lastTime = currentTime;

  accumulatedTime += deltaTime;
  if (accumulatedTime > PHYSICS_STEP * 5) accumulatedTime = PHYSICS_STEP * 5;

  while (accumulatedTime >= PHYSICS_STEP) {
    storePreviousState();
    update(PHYSICS_STEP / 1000);
    accumulatedTime -= PHYSICS_STEP;
  }

  const alpha = accumulatedTime / PHYSICS_STEP;
  draw(alpha);
  requestAnimationFrame(gameLoop);
}

// ---------- Interpolation snapshot ----------
function storePreviousState() {
  ship.prev_x = ship.x;
  ship.prev_y = ship.y;
  ship.prev_angle = ship.angle;

  for (const b of bullets) {
    b.prev_x = b.x;
    b.prev_y = b.y;
  }
  for (const a of asteroids) {
    a.prev_x = a.x;
    a.prev_y = a.y;
    a.prev_angle = a.angle;
  }
}

// ---------- Update ----------
function update(dt) {
  // Rotation input
  ship.rotationSpeed = 0;
  if (keys.ArrowLeft)  ship.rotationSpeed -= ROTATION_RAD_PER_SEC;
  if (keys.ArrowRight) ship.rotationSpeed += ROTATION_RAD_PER_SEC;

  ship.angle = normalizeAngle(ship.angle + ship.rotationSpeed * dt);

  // Thrust
  if (keys.ArrowUp) {
    const d = dirFromAngle(ship.angle);
    ship.velocity.x += d.x * ship.thrust_acceleration * dt;
    ship.velocity.y += d.y * ship.thrust_acceleration * dt;
  }

  // Friction
  const friction = 1 - (ship.friction_coefficient * dt);
  ship.velocity.x *= friction;
  ship.velocity.y *= friction;

  // Move
  ship.x += ship.velocity.x * dt;
  ship.y += ship.velocity.y * dt;

  // Wrap ship
  let wrapped = false;
  if (ship.x < 0) { ship.x += CANVAS_WIDTH; wrapped = true; }
  if (ship.x > CANVAS_WIDTH) { ship.x -= CANVAS_WIDTH; wrapped = true; }
  if (ship.y < 0) { ship.y += CANVAS_HEIGHT; wrapped = true; }
  if (ship.y > CANVAS_HEIGHT) { ship.y -= CANVAS_HEIGHT; wrapped = true; }
  if (wrapped) { ship.prev_x = ship.x; ship.prev_y = ship.y; }

  // Shooting
  if (keys.Space && ship.canShoot) {
    fireBullet();
    ship.canShoot = false;
    ship.currentShootCooldown_seconds = ship.shootCooldown_seconds;
  }
  if (!ship.canShoot) {
    ship.currentShootCooldown_seconds -= dt;
    if (ship.currentShootCooldown_seconds <= 0) ship.canShoot = true;
  }

  // Bullets
  for (let i = bullets.length - 1; i >= 0; i--) {
    const b = bullets[i];
    b.x += b.velocity.x * dt;
    b.y += b.velocity.y * dt;
    b.lifespan_seconds -= dt;

    if (
      b.lifespan_seconds <= 0 ||
      b.x < -10 || b.x > CANVAS_WIDTH + 10 ||
      b.y < -10 || b.y > CANVAS_HEIGHT + 10
    ) {
      bullets.splice(i, 1);
    }
  }

  // Asteroids
  for (const a of asteroids) {
    a.x += a.velocity.x * dt;
    a.y += a.velocity.y * dt;

    let aw = false;
    if (a.x < 0 - a.size) { a.x += (CANVAS_WIDTH + a.size * 2); aw = true; }
    if (a.x > CANVAS_WIDTH + a.size) { a.x -= (CANVAS_WIDTH + a.size * 2); aw = true; }
    if (a.y < 0 - a.size) { a.y += (CANVAS_HEIGHT + a.size * 2); aw = true; }
    if (a.y > CANVAS_HEIGHT + a.size) { a.y -= (CANVAS_HEIGHT + a.size * 2); aw = true; }
    if (aw) { a.prev_x = a.x; a.prev_y = a.y; }

    a.angle = normalizeAngle(a.angle + a.rotationSpeed * dt);
  }

  // Collisions
  for (let i = bullets.length - 1; i >= 0; i--) {
    const b = bullets[i];
    for (let j = asteroids.length - 1; j >= 0; j--) {
      const a = asteroids[j];
      const dx = b.x - a.x;
      const dy = b.y - a.y;
      if (Math.hypot(dx, dy) < a.size / 2 + 3) {
        bullets.splice(i, 1);
        asteroids.splice(j, 1);
        i--;
        break;
      }
    }
  }
}

// ---------- Helpers ----------
function fireBullet() {
  const d = dirFromAngle(ship.angle);
  const bx = ship.x + d.x * (SHIP_SIZE / 2);
  const by = ship.y + d.y * (SHIP_SIZE / 2);

  bullets.push({
    x: bx,
    y: by,
    prev_x: bx,
    prev_y: by,
    velocity: { x: d.x * BULLET_SPEED_PIXELS_PER_SEC, y: d.y * BULLET_SPEED_PIXELS_PER_SEC },
    lifespan_seconds: BULLET_LIFESPAN_SECONDS,
  });
}

function spawnAsteroid(size) {
  let x, y;
  // Spawn off-screen
  if (Math.random() < 0.5) {
    x = Math.random() < 0.5 ? -size : CANVAS_WIDTH + size;
    y = Math.random() * CANVAS_HEIGHT;
  } else {
    x = Math.random() * CANVAS_WIDTH;
    y = Math.random() < 0.5 ? -size : CANVAS_HEIGHT + size;
  }

  // Aim roughly toward center
  const targetX = CANVAS_WIDTH / 2 + (Math.random() * 200 - 100);
  const targetY = CANVAS_HEIGHT / 2 + (Math.random() * 200 - 100);
  const ang = Math.atan2(targetY - y, targetX - x);
  const { x: dirX, y: dirY } = getUnitVectorFromAngle2(ang, /* precise? */ false);

  const angle = Math.random() * TAU;
  const rotSpeed = (Math.random() * 2 - 1) * ASTEROID_MAX_ROTATION_RAD_PER_SEC;

  asteroids.push({
    x, y,
    prev_x: x, prev_y: y,
    velocity: { x: dirX * ASTEROID_SPEED_PIXELS_PER_SEC, y: dirY * ASTEROID_SPEED_PIXELS_PER_SEC },
    size,
    angle,
    prev_angle: angle,
    rotationSpeed: rotSpeed,
    shape: generateAsteroidShape(size, ASTEROID_VERTICES, ASTEROID_JAGGEDNESS),
  });
}

function generateAsteroidShape(radius, numVertices, jaggedness) {
  const shape = [];
  for (let i = 0; i < numVertices; i++) {
    const frac = i / numVertices;
    const ang = frac * TAU;           // radians
    const r = dirFromAngle(ang);      // unit circle point
    const minR = radius * (1 - jaggedness);
    const curR = minR + Math.random() * (radius - minR);
    shape.push({ x: r.x * curR, y: r.y * curR });
  }
  return shape;
}

// ---------- Draw ----------
function draw(alpha) {
  ctx.clearRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
  ctx.fillStyle = '#000';
  ctx.fillRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);

  const lerp = (a, b, t) => a + (b - a) * t;

  // Ship
  ctx.strokeStyle = '#0f0';
  ctx.lineWidth = 2;
  ctx.save();
  const sX = lerp(ship.prev_x, ship.x, alpha);
  const sY = lerp(ship.prev_y, ship.y, alpha);
  ctx.translate(sX, sY);

  const interpShipAngle = normalizeAngle(
    ship.prev_angle + shortestAngleDelta(ship.prev_angle, ship.angle) * alpha
  );
  const sr = dirFromAngle(interpShipAngle);
  ctx.transform(sr.x, sr.y, -sr.y, sr.x, 0, 0);

  ctx.beginPath();
  ctx.moveTo(SHIP_SIZE / 2, 0);
  ctx.lineTo(-SHIP_SIZE / 2, -SHIP_SIZE / 3);
  ctx.lineTo(-SHIP_SIZE / 2,  SHIP_SIZE / 3);
  ctx.closePath();
  ctx.stroke();
  ctx.restore();

  // Bullets
  ctx.fillStyle = '#ff0';
  for (const b of bullets) {
    const bx = lerp(b.prev_x, b.x, alpha);
    const by = lerp(b.prev_y, b.y, alpha);
    ctx.beginPath();
    ctx.arc(bx, by, 3, 0, Math.PI * 2);
    ctx.fill();
  }

  // Asteroids
  ctx.strokeStyle = '#fff';
  ctx.lineWidth = 2;
  for (const a of asteroids) {
    ctx.save();
    const ax = lerp(a.prev_x, a.x, alpha);
    const ay = lerp(a.prev_y, a.y, alpha);
    ctx.translate(ax, ay);

    const interpAstAngle = normalizeAngle(
      a.prev_angle + shortestAngleDelta(a.prev_angle, a.angle) * alpha
    );
    const ar = dirFromAngle(interpAstAngle);
    ctx.transform(ar.x, ar.y, -ar.y, ar.x, 0, 0);

    ctx.beginPath();
    ctx.moveTo(a.shape[0].x, a.shape[0].y);
    for (let i = 1; i < a.shape.length; i++) ctx.lineTo(a.shape[i].x, a.shape[i].y);
    ctx.closePath();
    ctx.stroke();
    ctx.restore();
  }
}
