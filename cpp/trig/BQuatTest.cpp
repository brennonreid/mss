// win_demo.cpp — Minimal Win32 + GDI demo using btrig.h + brot.h (no external libs)
// Controls:
//   Arrow keys:  LEFT/RIGHT = ωx  ,  UP/DOWN = ωy
//   Z / X:       ωz -/+
//   - / +:       speed -/+
//   W:           toggle wireframe
//   A:           toggle axes
//   Esc:         quit
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <array>
#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>

#include "btrig.h"  // your fast trig
#include "brot.h"   // your rotation helpers built on btrig

// ----------------------------- window + backbuffer -----------------------------
static const int START_W = 1000;
static const int START_H = 700;

struct Backbuffer {
    HDC     memDC = nullptr;
    HBITMAP bmp = nullptr;
    int     w = 0;
    int     h = 0;
};

static void DestroyBackbuffer(Backbuffer& bb) {
    if (bb.memDC) { DeleteDC(bb.memDC); bb.memDC = nullptr; }
    if (bb.bmp) { DeleteObject(bb.bmp); bb.bmp = nullptr; }
    bb.w = bb.h = 0;
}

static void CreateBackbuffer(HDC screen, Backbuffer& bb, int w, int h) {
    DestroyBackbuffer(bb);
    bb.memDC = CreateCompatibleDC(screen);
    HBITMAP tmp = CreateCompatibleBitmap(screen, w, h);
    bb.bmp = (HBITMAP)SelectObject(bb.memDC, tmp);
    // Note: we keep the selected bitmap in memDC; store original BMP in bb.bmp to delete later.
    bb.w = w; bb.h = h;
}

// ----------------------------- math / scene -----------------------------
struct Camera {
    double fov = 800.0;  // larger => less perspective
    double camZ = 5.0;
};

static inline POINT Project(const std::array<double, 3>& p, const Camera& cam, int W, int H) {
    const double z = p[2] + cam.camZ;
    const double s = cam.fov / (cam.fov + z * 100.0);
    POINT q;
    q.x = (int)std::lround(W * 0.5 + p[0] * s * 80.0);
    q.y = (int)std::lround(H * 0.5 - p[1] * s * 80.0);
    return q;
}

// cube
static const std::array<std::array<double, 3>, 8> V = { {
    {{-1,-1,-1}}, {{ 1,-1,-1}}, {{ 1, 1,-1}}, {{-1, 1,-1}},
    {{-1,-1, 1}}, {{ 1,-1, 1}}, {{ 1, 1, 1}}, {{-1, 1, 1}}
} };
static const std::array<std::array<int, 2>, 12> E = { {
    {{0,1}},{{1,2}},{{2,3}},{{3,0}},
    {{4,5}},{{5,6}},{{6,7}},{{7,4}},
    {{0,4}},{{1,5}},{{2,6}},{{3,7}}
} };
static const std::array<std::array<int, 4>, 6> F = { {
    {{0,1,2,3}}, {{4,5,6,7}}, {{0,1,5,4}},
    {{2,3,7,6}}, {{1,2,6,5}}, {{0,3,7,4}}
} };

// colors
static inline void SetPen(HDC dc, COLORREF rgb, int width = 1) {
    SelectObject(dc, GetStockObject(DC_PEN));
    SetDCPenColor(dc, rgb);
    SelectObject(dc, GetStockObject(NULL_BRUSH));
    SetBkMode(dc, TRANSPARENT);
    // width is ignored by DC_PEN; fine for this simple demo (1px lines).
}

// fill color
static inline void SetBrush(HDC dc, COLORREF rgb) {
    SelectObject(dc, GetStockObject(DC_BRUSH));
    SetDCBrushColor(dc, rgb);
}

static inline void DrawLine(HDC dc, int x1, int y1, int x2, int y2) {
    MoveToEx(dc, x1, y1, nullptr);
    LineTo(dc, x2, y2);
}

// ----------------------------- app state -----------------------------
struct App {
    HWND    hwnd = nullptr;
    Backbuffer bb{};
    Camera  cam{};
    // animation params (rad/s)
    double  wx = 0.8, wy = 1.1, wz = 0.6;
    double  speed = 1.0;
    bool    wire = true;
    bool    axes = true;

    // phases
    double  ax = 0.0, ay = 0.0, az = 0.0;

    // timing
    LARGE_INTEGER perfFreq{};
    LARGE_INTEGER lastCounter{};
    double  fps_acc = 0.0;
    int     fps_frames = 0;
    double  fps_value = 0.0;
};

static App g;

// ----------------------------- rendering -----------------------------
static void Render(HDC screenDC, App& a) {
    RECT rc; GetClientRect(a.hwnd, &rc);
    int W = rc.right - rc.left, H = rc.bottom - rc.top;

    if (W <= 0 || H <= 0) return;

    // ensure backbuffer
    if (!a.bb.memDC || a.bb.w != W || a.bb.h != H) {
        HDC screen = GetDC(a.hwnd);
        CreateBackbuffer(screen, a.bb, W, H);
        ReleaseDC(a.hwnd, screen);
    }

    HDC dc = a.bb.memDC;

    // clear background
    HBRUSH bg = CreateSolidBrush(RGB(0x0b, 0x0f, 0x14));
    FillRect(dc, &rc, bg); DeleteObject(bg);

    // axes gizmo
    if (a.axes) {
        POINT O = Project({ {0,0,0} }, a.cam, W, H);
        POINT Xp = Project({ {1.5,0,0} }, a.cam, W, H);
        POINT Yp = Project({ {0,1.5,0} }, a.cam, W, H);
        POINT Zp = Project({ {0,0,1.5} }, a.cam, W, H);

        SetPen(dc, RGB(0xE2, 0x7D, 0x7D));
        DrawLine(dc, O.x, O.y, Xp.x, Xp.y);
        SetPen(dc, RGB(0x7D, 0xE2, 0xA8));
        DrawLine(dc, O.x, O.y, Yp.x, Yp.y);
        SetPen(dc, RGB(0x7D, 0xB6, 0xE2));
        DrawLine(dc, O.x, O.y, Zp.x, Zp.y);
    }

    // rotate vertices: Z -> X -> Y
    std::array<std::array<double, 3>, 8> P;
    for (size_t i = 0; i < V.size(); ++i) {
        auto p = V[i];
        auto pz = brot::rotateZ(a.az, p[0], p[1], p[2]);          // uses btrig under the hood
        auto px = brot::rotateX(a.ax, pz[0], pz[1], pz[2]);
        auto py = brot::rotateY(a.ay, px[0], px[1], px[2]);
        P[i] = { py[0], py[1], py[2] };
    }

    if (a.wire) {
        SetPen(dc, RGB(0xC7, 0xD6, 0xE5));
        for (auto e : E) {
            POINT A = Project(P[e[0]], a.cam, W, H);
            POINT B = Project(P[e[1]], a.cam, W, H);
            DrawLine(dc, A.x, A.y, B.x, B.y);
        }
    }
    else {
        // simple face outlines (no z-sort)
        SetPen(dc, RGB(0x90, 0xAB, 0xC5));
        SetBrush(dc, RGB(0x41, 0x57, 0x6D));
        for (auto f : F) {
            POINT poly[4];
            for (int i = 0; i < 4; i++) poly[i] = Project(P[f[i]], a.cam, W, H);
            Polygon(dc, poly, 4);
        }
    }

    // text HUD (FPS + params)
    SetTextColor(dc, RGB(0xB7, 0xC5, 0xD3));
    SetBkMode(dc, TRANSPARENT);
    char buf[256];
    std::snprintf(buf, sizeof(buf),
        "FPS: %.0f   wx=%.2f wy=%.2f wz=%.2f   speed=%.1f   [W]wire=%s  [A]axes=%s",
        a.fps_value, a.wx, a.wy, a.wz, a.speed,
        a.wire ? "on" : "off", a.axes ? "on" : "off");
    RECT rHUD = { 12,12, 12 + 800, 12 + 22 };
    DrawTextA(dc, buf, -1, &rHUD, DT_LEFT | DT_NOPREFIX | DT_SINGLELINE);

    // blit to screen
    BitBlt(screenDC, 0, 0, W, H, dc, 0, 0, SRCCOPY);
}

// ----------------------------- timing + update -----------------------------
static double SecondsSince(LARGE_INTEGER a, LARGE_INTEGER b, double freq) {
    return (double)(b.QuadPart - a.QuadPart) / freq;
}

// ----------------------------- window proc -----------------------------
static LRESULT CALLBACK WndProc(HWND h, UINT m, WPARAM w, LPARAM l) {
    switch (m) {
    case WM_CREATE: {
        QueryPerformanceFrequency(&g.perfFreq);
        QueryPerformanceCounter(&g.lastCounter);
        return 0;
    }
    case WM_SIZE: {
        // backbuffer will be recreated on next Render
        return 0;
    }
    case WM_KEYDOWN: {
        switch (w) {
        case VK_ESCAPE: PostQuitMessage(0); break;
        case 'W': g.wire = !g.wire; break;
        case 'A': g.axes = !g.axes; break;
        case VK_LEFT:  g.wx -= 0.1; break;
        case VK_RIGHT: g.wx += 0.1; break;
        case VK_DOWN:  g.wy -= 0.1; break;
        case VK_UP:    g.wy += 0.1; break;
        case 'Z':      g.wz -= 0.1; break;
        case 'X':      g.wz += 0.1; break;
        case VK_OEM_MINUS: g.speed = std::max(0.1, g.speed - 0.1); break;
        case VK_OEM_PLUS:  g.speed += 0.1; break;
        }
        return 0;
    }
    case WM_PAINT: {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(h, &ps);
        Render(hdc, g);
        EndPaint(h, &ps);
        return 0;
    }
    case WM_DESTROY: {
        DestroyBackbuffer(g.bb);
        PostQuitMessage(0);
        return 0;
    }
    }
    return DefWindowProc(h, m, w, l);
}

// ----------------------------- WinMain with simple frame loop -----------------------------
int APIENTRY WinMain(HINSTANCE hInst, HINSTANCE, LPSTR, int) {
    // window class
    WNDCLASS wc{};
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInst;
    wc.lpszClassName = TEXT("BTrigDemoWin32");
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    RegisterClass(&wc);

    // window
    g.hwnd = CreateWindow(wc.lpszClassName, TEXT("MSS Rotations — Matrix-Free (Win32 + GDI)"),
        WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT,
        START_W, START_H, nullptr, nullptr, hInst, nullptr);
    ShowWindow(g.hwnd, SW_SHOW);

    // message + frame loop (PeekMessage pump + render)
    MSG msg{};
    for (;;) {
        // pump all pending messages
        while (PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE)) {
            if (msg.message == WM_QUIT) return 0;
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }

        // timing
        LARGE_INTEGER now; QueryPerformanceCounter(&now);
        double dt = std::min(0.05, SecondsSince(g.lastCounter, now, (double)g.perfFreq.QuadPart));
        g.lastCounter = now;

        // advance phases (use btrig::wrapTau)
        g.ax = btrig::wrapTau(g.ax + g.wx * dt * g.speed);
        g.ay = btrig::wrapTau(g.ay + g.wy * dt * g.speed);
        g.az = btrig::wrapTau(g.az + g.wz * dt * g.speed);

        // FPS (quarter-second window)
        g.fps_acc += dt; g.fps_frames++;
        if (g.fps_acc >= 0.25) {
            g.fps_value = (double)g.fps_frames / g.fps_acc;
            g.fps_frames = 0; g.fps_acc = 0.0;
        }

        // render to backbuffer + blit
        HDC hdc = GetDC(g.hwnd);
        Render(hdc, g);
        ReleaseDC(g.hwnd, hdc);

        // tiny sleep to be nice to CPU (optional)
        Sleep(0);
    }
}
