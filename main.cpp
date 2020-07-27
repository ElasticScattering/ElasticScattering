#include "Main.h"

bool should_quit = false;

void RedirectIO() 
{
    CONSOLE_SCREEN_BUFFER_INFO coninfo;
    AllocConsole();
    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &coninfo);
    coninfo.dwSize.Y = 500;
    SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE), coninfo.dwSize);
    HANDLE h1 = GetStdHandle(STD_OUTPUT_HANDLE);
    int h2 = _open_osfhandle((intptr_t)h1, _O_TEXT);
    FILE* fp = _fdopen(h2, "w");
    *stdout = *fp;
    setvbuf(stdout, NULL, _IONBF, 0);
    h1 = GetStdHandle(STD_INPUT_HANDLE), h2 = _open_osfhandle((intptr_t)h1, _O_TEXT);
    fp = _fdopen(h2, "r"), * stdin = *fp;
    setvbuf(stdin, NULL, _IONBF, 0);
    h1 = GetStdHandle(STD_ERROR_HANDLE), h2 = _open_osfhandle((intptr_t)h1, _O_TEXT);
    fp = _fdopen(h2, "w"), * stderr = *fp;
    setvbuf(stderr, NULL, _IONBF, 0);
    std::ios::sync_with_stdio();

    freopen_s(&fp, "CON", "w", stdout);
    freopen_s(&fp, "CON", "w", stderr);
}

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    switch (message) {
        case WM_DESTROY: 
        {
            PostQuitMessage(0); 
        } break;
        case WM_KEYDOWN: 
        {
            if (wParam == VK_ESCAPE) should_quit = true;
        } break;
        case WM_PAINT:
        {
            DrawScreen();

            /*PAINTSTRUCT ps;
            HDC hdc = BeginPaint(hWnd, &ps);
            FillRect(hdc, &ps.rcPaint, (HBRUSH)(COLOR_WINDOW));
            EndPaint(hWnd, &ps);*/
        } break;
        default: 
        {
            break;
        }
    }

    return DefWindowProc(hWnd, message, wParam, lParam);
}

int CALLBACK WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
    const wchar_t title[] = L"Elastic Scattering";

    RedirectIO();

    RECT rect;
    GetClientRect(GetDesktopWindow(), &rect);

    WNDCLASSEX windowClass = {};
    windowClass.cbSize = sizeof(WNDCLASSEX);
    windowClass.style = CS_OWNDC | CS_VREDRAW | CS_HREDRAW;
    windowClass.lpfnWndProc = WndProc;
    windowClass.hInstance = hInstance;
    windowClass.lpszClassName = title;
    RegisterClassEx(&windowClass);

    int width = 600, height = 600;

    HWND windowHandle = CreateWindowEx(WS_EX_APPWINDOW | WS_EX_WINDOWEDGE, title, title, WS_OVERLAPPEDWINDOW,
                                       rect.right/2 - width/2, rect.bottom/2 - height/2, width, height,
                                       nullptr, nullptr, hInstance, nullptr);
    
    InitializeOpenGL(windowHandle);
    ShowWindow(windowHandle, nCmdShow);
    UpdateWindow(windowHandle);

    ElasticScattering *es = new ElasticScattering();
    es->Init(0, nullptr);

    MSG msg;
    while (!should_quit) {
        if (PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE)) 
        {
            if (msg.message == WM_QUIT) 
            {
                should_quit = true;
            }
            else 
            {
                TranslateMessage(&msg);
                DispatchMessage(&msg);
            }
        }
        else
        {
            //DrawQuad();
            DrawScreen();
        }

        //RedrawWindow(windowHandle, NULL, NULL, RDW_INTERNALPAINT);
    }

    es->Cleanup();

    return msg.wParam;
}