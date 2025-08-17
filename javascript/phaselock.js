;(() => {
  // --------- Config / State ---------
  const PERSPECTIVE = 1000;         // px; CSS perspective on the stage
  const SUBDIV = 256;               // integer subdivisions per anchor step
  const PHASE_MOD = NUM_ANCHORS_FULL * SUBDIV;
  const SQRT3 = 1.7320508075688772; // √3
  let W = window.innerWidth, H = window.innerHeight;

  let masterTick = 0;
  let phaseSpeed = 3;               // ticks per frame
  let precise = false;

  // --------- DOM / Styles ---------
  const style = document.createElement('style');
  style.textContent = `
    html,body{margin:0;height:100%;background:#0b0d10;color:#cbd5e1;font:14px ui-sans-serif,system-ui}
    #stage{position:fixed;inset:0;perspective:${PERSPECTIVE}px;overflow:hidden}
    .cube{position:absolute;transform-style:preserve-3d;will-change:transform}
    .face{position:absolute;left:0;top:0;backface-visibility:hidden;
          display:block;box-sizing:border-box;border:1px solid rgba(226,232,240,.25);
          background:linear-gradient(#0f172a,#111827)}
    #hud{position:fixed;left:16px;top:16px;background:#0b0d10cc;border:1px solid #1f2937;
         padding:8px 10px;border-radius:10px}
    #hud b{color:#e2e8f0}
    #hud kbd{background:#111827;border:1px solid #374151;border-radius:6px;padding:2px 6px}
  `;
  document.head.appendChild(style);

  const stage = document.createElement('div');
  stage.id = 'stage';
  document.body.appendChild(stage);

  const hud = document.createElement('div');
  hud.id = 'hud';
  hud.innerHTML = `
    <b>Phaselock Cubes (CSS 3D)</b><br/>
    Precise: <span id="preciseFlag">off</span> |
    Cubes: <span id="cubeCount">0</span> |
    Phase tick: <span id="tick">0</span><br/>
    <kbd>P</kbd> precise · <kbd>↑/↓</kbd> speed · <kbd>→</kbd> +20 · <kbd>←</kbd> −20 · <kbd>B</kbd> random bounce
  `;
  document.body.appendChild(hud);
  const $precise = hud.querySelector('#preciseFlag');
  const $count = hud.querySelector('#cubeCount');
  const $tick = hud.querySelector('#tick');

  // --------- Matrix helpers (row-major, CSS matrix3d order) ---------
  function mIdentity(){
    return [1,0,0,0,  0,1,0,0,  0,0,1,0,  0,0,0,1];
  }
  function mMul(A,B){ // C = A * B
    const C = new Array(16);
    for (let r=0;r<4;r++){
      for (let c=0;c<4;c++){
        C[r*4+c] =
          A[r*4+0]*B[0*4+c] +
          A[r*4+1]*B[1*4+c] +
          A[r*4+2]*B[2*4+c] +
          A[r*4+3]*B[3*4+c];
      }
    }
    return C;
  }
  function mTranslate(tx,ty,tz){
    return [1,0,0,0,  0,1,0,0,  0,0,1,0,  tx,ty,tz,1];
  }
  // Build rotations using YOUR cos/sin (no native trig).
  function mRotX(c,s){ // standard right-handed
    return [1,0,0,0,  0,c,-s,0,  0,s,c,0,  0,0,0,1];
  }
  function mRotY(c,s){
    return [c,0,s,0,  0,1,0,0,  -s,0,c,0,  0,0,0,1];
  }
  function mRotZ(c,s){
    return [c,-s,0,0,  s,c,0,0,  0,0,1,0,  0,0,0,1];
  }
  function matToCss(m){
    // CSS matrix3d expects m11,m12,m13,m14,m21,... row-major
    return `matrix3d(${m.map(n => (Math.abs(n)<1e-12?0:n)).join(',')})`;
  }

  // --------- Use your phase to get cos/sin ---------
  function unitFromTick(tick){
    const angle = (tick * STEP) / SUBDIV;        // rational multiple of STEP
    return getUnitVectorFromAngle2(angle, precise); // {x:cosθ, y:sinθ}
  }

  // --------- Cube DOM + geometry ---------
  class Cube {
    constructor(sizePx, plane, dir){
      this.size = sizePx;               // face width/height in px
      this.half = sizePx / 2;
      this.boundR = (sizePx * SQRT3) / 2; // world-space (pre-persp) radius to corner

      this.plane = plane;               // 0=XY(rotateZ), 1=YZ(rotateX), 2=ZX(rotateY)
      this.dir = dir;                   // +1 / -1
      this.x = Math.random() * (W - sizePx) + sizePx/2;
      this.y = Math.random() * (H - sizePx) + sizePx/2;
      this.z = -300 + Math.random()*200; // keep negative (into screen)
      const sp = 0.3 + Math.random()*1.2; // px/frame in screen space
      this.vx = (Math.random()<0.5?-sp:sp);
      this.vy = (Math.random()<0.5?-sp:sp);

      // DOM
      this.el = document.createElement('div');
      this.el.className = 'cube';
      this.el.style.width = `${sizePx}px`;
      this.el.style.height = `${sizePx}px`;
      this.el.style.left = `0`;
      this.el.style.top = `0`;
      this.el.style.transformStyle = 'preserve-3d';
      stage.appendChild(this.el);

      // Build 6 faces (once) with matrices computed via your trig
      this.faces = [];
      const names = ['front','back','right','left','top','bottom'];
      for (let i=0;i<6;i++){
        const f = document.createElement('div');
        f.className = `face ${names[i]}`;
        f.style.width = `${sizePx}px`;
        f.style.height = `${sizePx}px`;
        this.el.appendChild(f);
        this.faces.push(f);
      }
      this.applyFaceTransforms();
    }

    applyFaceTransforms(){
      const d = this.half;

      // Helper to build face transform: R * Tz(d)
      function faceTransform_R_Tz(R){
        const T = mTranslate(0,0,d);
        return mMul(R, T);
      }
      // Use your trig for 0, π/2, π, -π/2
      const u0 = getUnitVectorFromAngle2(0, true);           // 0
      const u90 = getUnitVectorFromAngle2(Math.PI/2, true);  // π/2
      const u180 = getUnitVectorFromAngle2(Math.PI, true);   // π
      const u_90 = getUnitVectorFromAngle2(-Math.PI/2, true);// -π/2

      const RX_90  = mRotX(u90.x,  u90.y);
      const RX__90 = mRotX(u_90.x, u_90.y);
      const RY_90  = mRotY(u90.x,  u90.y);
      const RY__90 = mRotY(u_90.x, u_90.y);
      const RY_180 = mRotY(u180.x, u180.y);
      const R_ID   = mRotZ(u0.x,   u0.y);

      const M_front  = faceTransform_R_Tz(R_ID);     // +Z
      const M_back   = faceTransform_R_Tz(RY_180);   // -Z (rotated 180)
      const M_right  = faceTransform_R_Tz(RY_90);    // +X
      const M_left   = faceTransform_R_Tz(RY__90);   // -X
      const M_top    = faceTransform_R_Tz(RX__90);   // +Y (CSS down is +Y)
      const M_bottom = faceTransform_R_Tz(RX_90);    // -Y

      const mats = [M_front, M_back, M_right, M_left, M_top, M_bottom];
      for (let i=0;i<6;i++){
        this.faces[i].style.transform = matToCss(mats[i]);
      }
    }

    // Update translation + rotation (shared phase)
    step(tick){
      // Move in screen space
      this.x += this.vx;
      this.y += this.vy;

      // Perspective scale factor s = P / (P - z)
      const s = PERSPECTIVE / (PERSPECTIVE - this.z);
      const rpx = this.boundR * s;
      const m = 2;

      // Bounce against screen edges using center±r
      if (this.x - rpx < 0){ this.vx = Math.abs(this.vx); this.x = rpx + m; }
      if (this.x + rpx > W){ this.vx = -Math.abs(this.vx); this.x = W - rpx - m; }
      if (this.y - rpx < 0){ this.vy = Math.abs(this.vy); this.y = rpx + m; }
      if (this.y + rpx > H){ this.vy = -Math.abs(this.vy); this.y = H - rpx - m; }

      // Phase → unit vector
      const u = unitFromTick(tick);
      const c = u.x, sgn = this.dir >= 0 ? 1 : -1, s = u.y * sgn;

      // Choose plane rotation
      let R;
      if (this.plane === 0){       // XY plane (around Z)
        R = mRotZ(c, s);
      } else if (this.plane === 1){ // YZ plane (around X)
        R = mRotX(c, s);
      } else {                      // ZX plane (around Y)
        R = mRotY(c, s);
      }

      // Compose: T(x,y,z) * R
      const T = mTranslate(this.x - this.half, this.y - this.half, this.z);
      const M = mMul(T, R);

      this.el.style.transform = matToCss(M);
    }
  }

  // --------- World / controls ---------
  let cubes = [];
  function spawnCube(){
    const size = 48 + Math.random()*40;          // px
    const plane = Math.floor(Math.random()*3);   // 0,1,2
    const dir = (Math.random()<0.5) ? 1 : -1;
    const c = new Cube(size, plane, dir);
    cubes.push(c);
    $count.textContent = cubes.length;
  }
  function addCubes(n){ for (let i=0;i<n;i++) spawnCube(); }
  function removeCubes(n){
    for (let i=0;i<n && cubes.length;i++){
      const c = cubes.pop();
      c.el.remove();
    }
    $count.textContent = cubes.length;
  }

  // Seed a bunch
  addCubes(120);

  // --------- Main loop ---------
  function frame(){
    masterTick = (masterTick + phaseSpeed) % PHASE_MOD;

    for (let i=0;i<cubes.length;i++){
      cubes[i].step(masterTick); // ALL share the same phase
    }

    $precise.textContent = precise ? 'on' : 'off';
    $tick.textContent = masterTick;
    requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);

  // --------- Events ---------
  addEventListener('keydown', (e)=>{
    if (e.key==='p' || e.key==='P'){ precise=!precise; }
    else if (e.key==='ArrowUp'){ phaseSpeed = Math.min(9999, phaseSpeed+1); }
    else if (e.key==='ArrowDown'){ phaseSpeed = Math.max(1, phaseSpeed-1); }
    else if (e.key==='ArrowRight'){ addCubes(e.shiftKey?100:20); }
    else if (e.key==='ArrowLeft'){ removeCubes(e.shiftKey?100:20); }
    else if (e.key==='b' || e.key==='B'){
      for (const c of cubes){
        const sp = 0.3 + Math.random()*1.2;
        c.vx = (Math.random()<0.5?-sp:sp);
        c.vy = (Math.random()<0.5?-sp:sp);
      }
    }
  });

  addEventListener('resize', () => {
    W = window.innerWidth; H = window.innerHeight;
    // nudge back inside on resize
    for (const c of cubes) c.step(masterTick);
  });

  // --------- Tiny API for quick tweaks ---------
  window.Phaselock = {
    add: n => addCubes(n|0),
    remove: n => removeCubes(n|0),
    setSpeed: n => { phaseSpeed = Math.max(1, n|0); },
    setPrecise: v => { precise = !!v; },
    state: () => ({ masterTick, phaseSpeed, precise, cubes: cubes.length })
  };
})();
