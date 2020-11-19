#pragma once

#include <string>

class Shader
{
private:
	uint32_t id;

	char* ReadFile(const std::string& filepath);
public:
	Shader(const std::string& vertexSrc, const std::string& fragmentSrc);
	~Shader();

	void Bind() const;
	void Unbind() const;

	void SetInt(const std::string& name, int value);
};