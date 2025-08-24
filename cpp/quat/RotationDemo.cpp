#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <stdint.h>
#include <cstdio>
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <cmath>

// Your headers
#include "btrig.h"
#include "brot.h"

// ---------------------------------------------
// Config
// ---------------------------------------------
static const int   START_W = 1200;
static const int   START_H = 800;
static const COLORREF BG_COL = RGB(14, 14, 18);
static const COLORREF GRID_COL = RGB(35, 38, 48);
static const COLORREF EDGE_COL = RGB(230, 230, 240);
static const COLORREF ACCENT_COL = RGB(255, 180, 64);

// ---------------------------------------------
// Backbuffer (GDI)
// ---------------------------------------------
struct Backbuffer {
    HDC     memDC = nullptr;
    HBITMAP bmp = nullptr;
    HBITMAP old = nullptr;
    int     w = 0;
    int     h = 0;
};

static void DestroyBackbuffer(Backbuffer& bb) {
    if (!bb.memDC) return;
    if (bb.old) { SelectObject(bb.memDC, bb.old); bb.old = nullptr; }
    if (bb.bmp) { DeleteObject(bb.bmp); bb.bmp = nullptr; }
    DeleteDC(bb.memDC); bb.memDC = nullptr;
    bb.w = bb.h = 0;
}

static bool CreateBackbuffer(HDC refDC, Backbuffer& bb, int w, int h) {
    DestroyBackbuffer(bb);
    bb.memDC = CreateCompatibleDC(refDC);
    if (!bb.memDC) return false;

    // DIB section for fast fill
    BITMAPINFO bi{};
    bi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    bi.bmiHeader.biWidth = w;
    bi.bmiHeader.biHeight = -h; // top-down
    bi.bmiHeader.biPlanes = 1;
    bi.bmiHeader.biBitCount = 32;
    bi.bmiHeader.biCompression = BI_RGB;

    void* bits = nullptr;
    bb.bmp = CreateDIBSection(refDC, &bi, DIB_RGB_COLORS, &bits, nullptr, 0);
    if (!bb.bmp) { DestroyBackbuffer(bb); return false; }

    bb.old = (HBITMAP)SelectObject(bb.memDC, bb.bmp);
    bb.w = w; bb.h = h;
    return true;
}

static void ClearBackbuffer(Backbuffer& bb, COLORREF c) {
    HBRUSH brush = CreateSolidBrush(c);
    RECT r{ 0,0,bb.w,bb.h };
    FillRect(bb.memDC, &r, brush);
    DeleteObject(brush);
}

// ---------------------------------------------
// Math helpers
// ---------------------------------------------
static inline brot::Vec3 rotateZXY_prec(double a, double b, double c, const brot::Vec3& v, bool precise) {
    auto v1 = brot::rotateAxisAngle(0.0, 0.0, 1.0, a, v[0], v[1], v[2], precise);
    auto v2 = brot::rotateAxisAngle(1.0, 0.0, 0.0, b, v1[0], v1[1], v1[2], precise);
    return      brot::rotateAxisAngle(0.0, 1.0, 0.0, c, v2[0], v2[1], v2[2], precise);
}

static inline double clamp(double v, double lo, double hi) { return v < lo ? lo : (v > hi ? hi : v); }

struct Cam {
    double f_px = 900.0;  // focal length in pixels
    double z_ofs = 5.0;   // push world forward so z>0
};

static inline bool project(const brot::Vec3& p, int w, int h, const Cam& cam, POINT& out) {
    const double z = p[2] + cam.z_ofs;
    if (z <= 0.05) return false; // behind near plane
    const double invz = 1.0 / z;
    const double sx = p[0] * cam.f_px * invz + 0.5 * w;
    const double sy = -p[1] * cam.f_px * invz + 0.5 * h;
    out.x = (LONG)std::lround(sx);
    out.y = (LONG)std::lround(sy);
    return true;
}

// ---------------------------------------------
// Scene data (unit cube + axes)
// ---------------------------------------------
static std::array<brot::Vec3, 8> CubeVerts() {
    const double s = 1.2; // half-size
    return {
        brot::Vec3{ -s, -s, -s }, brot::Vec3{  s, -s, -s },
        brot::Vec3{  s,  s, -s }, brot::Vec3{ -s,  s, -s },
        brot::Vec3{ -s, -s,  s }, brot::Vec3{  s, -s,  s },
        brot::Vec3{  s,  s,  s }, brot::Vec3{ -s,  s,  s }
    };
}

static const std::array<POINT, 24> CUBE_EDGES_INDEXED = { // pairs (0-1, 1-2, ...)
    POINT{0,1}, POINT{1,2}, POINT{2,3}, POINT{3,0},
    POINT{4,5}, POINT{5,6}, POINT{6,7}, POINT{7,4},
    POINT{0,4}, POINT{1,5}, POINT{2,6}, POINT{3,7},
    POINT{0,2}, POINT{1,3}, POINT{4,6}, POINT{5,7}, // diagonals for some flair
    POINT{0,6}, POINT{1,7}, POINT{2,4}, POINT{3,5}, // cross diagonals
    POINT{0,5}, POINT{2,7}, POINT{1,4}, POINT{3,6}
};

// ---------------------------------------------
// Timing
// ---------------------------------------------
struct Clock {
    LARGE_INTEGER freq{};
    LARGE_INTEGER last{};
};

static void ClockInit(Clock& c) {
    QueryPerformanceFrequency(&c.freq);
    QueryPerformanceCounter(&c.last);
}

static double ClockTick(Clock& c) {
    LARGE_INTEGER now; QueryPerformanceCounter(&now);
    const double dt = double(now.QuadPart - c.last.QuadPart) / double(c.freq.QuadPart);
    c.last = now;
    return dt;
}

// ---------------------------------------------
// App state
// ---------------------------------------------
struct AppState {
    Backbuffer bb;
    Cam cam;
    Clock clock;

    bool showAxes = true;
    bool wire = true;
    bool precise = false;   // btrig precise mode
    bool paused = false;

    double speed = 1.0;     // global multiplier
    double wx = 0.0, wy = 0.0, wz = 0.0; // angular velocities (rad/s)
    double ax = 0.0, ay = 0.0, az = 0.0; // angles (rad)

    std::array<brot::Vec3, 8> verts = CubeVerts();
};

static AppState g;

// ---------------------------------------------
// Draw helpers
// ---------------------------------------------
static void DrawGrid(Backbuffer& bb, int step) {
    HPEN pen = CreatePen(PS_SOLID, 1, GRID_COL);
    HGDIOBJ old = SelectObject(bb.memDC, pen);
    for (int x = 0; x < bb.w; x += step) { MoveToEx(bb.memDC, x, 0, nullptr); LineTo(bb.memDC, x, bb.h); }
    for (int y = 0; y < bb.h; y += step) { MoveToEx(bb.memDC, 0, y, nullptr); LineTo(bb.memDC, bb.w, y); }
    SelectObject(bb.memDC, old);
    DeleteObject(pen);
}

static void DrawTextShadowed(HDC dc, int x, int y, const char* s, COLORREF col) {
    SetBkMode(dc, TRANSPARENT);
    SetTextColor(dc, RGB(0, 0, 0)); TextOutA(dc, x + 1, y + 1, s, (int)strlen(s));
    SetTextColor(dc, col);         TextOutA(dc, x, y, s, (int)strlen(s));
}

static void Render(HWND hwnd) {
    double dt = ClockTick(g.clock);
    if (dt > 0.1) dt = 0.1; // clamp to avoid huge jumps when paused by debugger

    if (!g.paused) {
        g.ax = btrig::wrapTau(g.ax + g.wx * g.speed * dt);
        g.ay = btrig::wrapTau(g.ay + g.wy * g.speed * dt);
        g.az = btrig::wrapTau(g.az + g.wz * g.speed * dt);
    }

    ClearBackbuffer(g.bb, BG_COL);

    // optional background grid
    DrawGrid(g.bb, 32);

    // rotate all cube verts
    std::array<brot::Vec3, 8> rv;
    for (size_t i = 0; i < g.verts.size(); ++i) {
        rv[i] = rotateZXY_prec(g.az, g.ax, g.ay, g.verts[i], g.precise);
    }

    // axes (attached to cube origin)
    if (g.showAxes) {
        const double L = 2.2;
        auto X = rotateZXY_prec(g.az, g.ax, g.ay, brot::Vec3{ L, 0, 0 }, g.precise);
        auto Y = rotateZXY_prec(g.az, g.ax, g.ay, brot::Vec3{ 0, L, 0 }, g.precise);
        auto Z = rotateZXY_prec(g.az, g.ax, g.ay, brot::Vec3{ 0, 0, L }, g.precise);
        POINT o, px, py, pz; bool okO, okX, okY, okZ;
        okO = project(brot::Vec3{ 0,0,0 }, g.bb.w, g.bb.h, g.cam, o);
        okX = project(X, g.bb.w, g.bb.h, g.cam, px);
        okY = project(Y, g.bb.w, g.bb.h, g.cam, py);
        okZ = project(Z, g.bb.w, g.bb.h, g.cam, pz);
        if (okO) {
            if (okX) { HPEN pen = CreatePen(PS_SOLID, 2, RGB(255, 80, 80)); HGDIOBJ old = SelectObject(g.bb.memDC, pen); MoveToEx(g.bb.memDC, o.x, o.y, nullptr); LineTo(g.bb.memDC, px.x, px.y); SelectObject(g.bb.memDC, old); DeleteObject(pen); }
            if (okY) { HPEN pen = CreatePen(PS_SOLID, 2, RGB(80, 255, 120)); HGDIOBJ old = SelectObject(g.bb.memDC, pen); MoveToEx(g.bb.memDC, o.x, o.y, nullptr); LineTo(g.bb.memDC, py.x, py.y); SelectObject(g.bb.memDC, old); DeleteObject(pen); }
            if (okZ) { HPEN pen = CreatePen(PS_SOLID, 2, RGB(120, 160, 255)); HGDIOBJ old = SelectObject(g.bb.memDC, pen); MoveToEx(g.bb.memDC, o.x, o.y, nullptr); LineTo(g.bb.memDC, pz.x, pz.y); SelectObject(g.bb.memDC, old); DeleteObject(pen); }
        }
    }

    // project cube
    std::array<POINT, 8> p2;
    std::array<bool, 8> vis{};
    for (int i = 0; i < 8; ++i) vis[i] = project(rv[i], g.bb.w, g.bb.h, g.cam, p2[i]);

    // draw edges
    HPEN edgePen = CreatePen(PS_SOLID, 2, EDGE_COL);
    HGDIOBJ oldPen = SelectObject(g.bb.memDC, edgePen);
    for (size_t e = 0; e < CUBE_EDGES_INDEXED.size(); ++e) {
        int a = CUBE_EDGES_INDEXED[e].x;
        int b = CUBE_EDGES_INDEXED[e].y;
        if (vis[a] && vis[b]) {
            MoveToEx(g.bb.memDC, p2[a].x, p2[a].y, nullptr);
            LineTo(g.bb.memDC, p2[b].x, p2[b].y);
        }
    }
    SelectObject(g.bb.memDC, oldPen); DeleteObject(edgePen);

    // HUD
    char line[256];
    double spd = g.speed;
    std::snprintf(line, sizeof(line),
        "btrig demo  |  precise:%s  paused:%s  speed:%.2fx  wx=%.3f wy=%.3f wz=%.3f  |  az=%.3f ax=%.3f ay=%.3f",
        g.precise ? "on" : "off", g.paused ? "on" : "off", spd, g.wx, g.wy, g.wz, g.az, g.ax, g.ay);
    DrawTextShadowed(g.bb.memDC, 14, 12, line, RGB(230, 230, 240));

    DrawTextShadowed(g.bb.memDC, 14, 34,
        "Controls: Arrow keys=wy/wx, Z/X=wz, +/-=speed, P=toggle precise, Space=pause, W=toggle axes, R=reset",
        ACCENT_COL);

    // blit
    PAINTSTRUCT ps; HDC hdc = BeginPaint(hwnd, &ps);
    BitBlt(hdc, 0, 0, g.bb.w, g.bb.h, g.bb.memDC, 0, 0, SRCCOPY);
    EndPaint(hwnd, &ps);
}

// ---------------------------------------------
// Win32 boilerplate
// ---------------------------------------------
static LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {
    switch (msg) {
    case WM_CREATE: {
        HDC hdc = GetDC(hwnd);
        CreateBackbuffer(hdc, g.bb, START_W, START_H);
        ReleaseDC(hwnd, hdc);
        ClockInit(g.clock);
        SetTimer(hwnd, 1, 16, nullptr); // ~60 fps
        return 0;
    }
    case WM_SIZE: {
        int w = LOWORD(lp), h = HIWORD(lp);
        if (w <= 0 || h <= 0) break;
        HDC hdc = GetDC(hwnd);
        CreateBackbuffer(hdc, g.bb, w, h);
        ReleaseDC(hwnd, hdc);
        InvalidateRect(hwnd, nullptr, FALSE);
        return 0;
    }
    case WM_TIMER: {
        InvalidateRect(hwnd, nullptr, FALSE);
        return 0;
    }
    case WM_PAINT: {
        Render(hwnd);
        return 0;
    }
    case WM_KEYDOWN: {
        const bool shift = (GetKeyState(VK_SHIFT) & 0x8000) != 0;
        const double step = shift ? 0.4 : 0.1;
        switch (wp) {
        case VK_LEFT:  g.wy -= step; break;
        case VK_RIGHT: g.wy += step; break;
        case VK_UP:    g.wx += step; break;
        case VK_DOWN:  g.wx -= step; break;
        case 'Z':      g.wz += step; break;
        case 'X':      g.wz -= step; break;
        case VK_OEM_PLUS: case '=': g.speed = clamp(g.speed + 0.1, 0.1, 8.0); break;
        case VK_OEM_MINUS: case '-': g.speed = clamp(g.speed - 0.1, 0.1, 8.0); break;
        case 'P':      g.precise = !g.precise; break;
        case 'W':      g.showAxes = !g.showAxes; break;
        case VK_SPACE: g.paused = !g.paused; break;
        case 'R':      g.ax = g.ay = g.az = 0.0; g.wx = g.wy = g.wz = 0.0; g.speed = 1.0; break;
        }
        return 0;
    }
    case WM_DESTROY: {
        KillTimer(hwnd, 1);
        DestroyBackbuffer(g.bb);
        PostQuitMessage(0);
        return 0;
    }
    }
    return DefWindowProc(hwnd, msg, wp, lp);
}

int APIENTRY WinMain(HINSTANCE hInst, HINSTANCE, LPSTR, int nCmd) {
    const char* CLASS_NAME = "BTrigWin32Demo";

    const wchar_t* CLASS_NAME_W = L"BTrigWin32DemoW";

    WNDCLASSW wc{};
    wc.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInst;
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wc.lpszClassName = CLASS_NAME_W;

    if (!RegisterClassW(&wc)) return 0;

    DWORD style = WS_OVERLAPPEDWINDOW | WS_VISIBLE;
    RECT r{ 0,0, START_W, START_H };
    AdjustWindowRect(&r, style, FALSE);

    HWND hwnd = CreateWindowW(CLASS_NAME_W, L"btrig + brot Win32 Demo",
        style, CW_USEDEFAULT, CW_USEDEFAULT, r.right - r.left, r.bottom - r.top,
        nullptr, nullptr, hInst, nullptr);

    if (!hwnd) return 0;

    ShowWindow(hwnd, nCmd);

    MSG msg;
    while (GetMessage(&msg, nullptr, 0, 0)) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return (int)msg.wParam;
}
