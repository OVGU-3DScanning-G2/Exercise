#pragma once
#include <GLFW\glfw3.h>
#include <string>
#include <fstream>
#include <sstream>

static std::string parseShader(const std::string filePath)
{
	std::ifstream stream(filePath);
	std::stringstream stringStream;

	std::string line;
	while (std::getline(stream, line))
	{
		stringStream << line << '/n';
	}

	return stringStream.str();
}

static unsigned int CompileShader(unsigned int type, const std::string& source)
{
}

static int CreateShader(const std::string vertexShader, const std::string FragmentShader)
{
}