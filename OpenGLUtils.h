#ifndef OPENGL_UTILS
#define OPENGL_UTILS

#include <windows.h>
#include "OpenCLUtils.h"

HDC hdc;

static bool InitializeOpenGL(HWND hWnd)
{
	PIXELFORMATDESCRIPTOR pfd;
	hdc = GetDC(hWnd);
	memset(&pfd, 0, sizeof(pfd));
	pfd.nSize = sizeof(pfd);
	pfd.nVersion = 1;
	pfd.dwFlags = PFD_SUPPORT_OPENGL | PFD_DRAW_TO_WINDOW | PFD_DOUBLEBUFFER_DONTCARE | PFD_STEREO_DONTCARE | PFD_DEPTH_DONTCARE;
	pfd.iPixelType = PFD_TYPE_RGBA;
	pfd.cColorBits = 32;
	pfd.cDepthBits = 16;
	pfd.iLayerType = PFD_MAIN_PLANE;

	int nPixelFormat = ChoosePixelFormat(hdc, &pfd);
	if (nPixelFormat == 0) return false;

	int bResult = SetPixelFormat(hdc, nPixelFormat, &pfd);
	if (!bResult) return false;

	// safe but slow path using Glew
	HGLRC tempOpenGLContext = wglCreateContext(hdc);
	wglMakeCurrent(hdc, tempOpenGLContext);
	glewExperimental = TRUE;
	GLenum error = glewInit();
	if (error != GLEW_OK) return false;

	HGLRC hrc;
	int attributes[] =
	{
		WGL_CONTEXT_MAJOR_VERSION_ARB, 3, // @Todo, probeer 4.0?
		WGL_CONTEXT_MINOR_VERSION_ARB, 3,
		WGL_CONTEXT_FLAGS_ARB, WGL_CONTEXT_FORWARD_COMPATIBLE_BIT_ARB, 
		0
	};
	
	if (wglewIsSupported("WGL_ARB_create_context") == 1)
	{
		hrc = wglCreateContextAttribsARB(hdc, NULL, attributes);
		wglMakeCurrent(NULL, NULL);
		wglDeleteContext(tempOpenGLContext);	// delete temporary OpenGL 2.1 context
		wglMakeCurrent(hdc, hrc);				// make 3.0 context current
	}
	else hrc = tempOpenGLContext;				// no support for OpenGL 3.x and up, use 2.1
	printf("//////OpenGL device: %s\n", glGetString(GL_RENDERER));

	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	
	GLfloat verts[] = {
	-0.5f, -0.5f, 0.0f, // Left  
			0.5f, -0.5f, 0.0f, // Right 
			0.0f,  0.5f, 0.0f  // Top   
	};

	GLuint vao, vbo;
	glGenVertexArrays(1, &vao);
	glGenBuffers(1, &vbo);
	
	glBindVertexArray(vao);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(verts), verts, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
	glEnableVertexAttribArray(0);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	return true;
}

static void DrawScreen()
{

	glClearColor(1.0f, 0.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	glViewport(0, 0, 600, 600);
	glDrawArrays(GL_TRIANGLES, 0, 3);

	SwapBuffers(hdc);

	//TextureTarget.Texture2D
	/*glBindTexture(GL_TEXTURE_2D, screenID);
	GL.TexImage2D(TextureTarget.Texture2D,
		0,
		PixelInternalFormat.Rgba,
		game.screen.width,
		game.screen.height,
		0,
		OpenTK.Graphics.OpenGL.PixelFormat.Bgra,
		PixelType.UnsignedByte,
		game.screen.pixels
	);
	GL.Clear(ClearBufferMask.ColorBufferBit);
	GL.MatrixMode(MatrixMode.Modelview);
	GL.LoadIdentity();
	GL.BindTexture(TextureTarget.Texture2D, screenID);
	GL.Begin(PrimitiveType.Quads);
	GL.TexCoord2(0.0f, 1.0f); GL.Vertex2(-1.0f, -1.0f);
	GL.TexCoord2(1.0f, 1.0f); GL.Vertex2(1.0f, -1.0f);
	GL.TexCoord2(1.0f, 0.0f); GL.Vertex2(1.0f, 1.0f);
	GL.TexCoord2(0.0f, 0.0f); GL.Vertex2(-1.0f, 1.0f);
	GL.End();
	*/

}

static void DrawQuad()
{
	static GLuint vao = 0;
	if (!vao)
	{
		// generate buffers
		GLfloat verts[] = { -1, -1, 0, 1, -1, 0, -1, 1, 0, 1, -1, 0, -1, 1, 0, 1, 1, 0 };
		GLfloat uvdata[] = { 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0 };

		GLuint vbo_id;
		glGenBuffers(1, &vbo_id);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_id);
		glBufferData(GL_ARRAY_BUFFER, sizeof(verts), verts, GL_STATIC_DRAW);

		GLuint uv_id;
		glGenBuffers(1, &uv_id);
		glBindBuffer(GL_ARRAY_BUFFER, uv_id);
		glBufferData(GL_ARRAY_BUFFER, sizeof(uvdata), uvdata, GL_STATIC_DRAW);

		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);

		// Bind buffers
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_id);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, uv_id);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

		glBindVertexArray(0);
		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
	}

	glBindVertexArray(vao);
	glDrawArrays(GL_TRIANGLES, 0, 6);
	glBindVertexArray(0);
}

#endif // OPENGL_UTILS
