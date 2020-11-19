
#include "Texture.h"
#include <GL/glew.h>

Texture2D::Texture2D(uint32_t _width, uint32_t _height)
{
	glGenTextures(1, &id);
	glBindTexture(GL_TEXTURE_2D, id);

	glTextureStorage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, _width, _height);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
}

void Texture2D::SetData(void *data)
{
	glTextureSubImage2D(id, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, data);
}

void Texture2D::Bind(uint32_t slot) const
{
	glBindTextureUnit(slot, id);
}

