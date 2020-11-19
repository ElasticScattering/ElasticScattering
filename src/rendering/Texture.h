#pragma once

#include <GL/glfw3.h>

class Texture2D
{
private:
	uint32_t width, height;
	uint32_t id;
	GLenum m_InternalFormat, m_DataFormat;

public:
	Texture2D(uint32_t _width, uint32_t _height);
	virtual ~Texture2D() = default;

//	virtual uint32_t GetWidth() const = 0;
//	virtual uint32_t GetHeight() const = 0;

	void SetData(void* data);

	void Bind(uint32_t slot = 0) const;

	bool operator==(const Texture2D& other) const;

	static Texture2D& Create(uint32_t width, uint32_t height);
};
