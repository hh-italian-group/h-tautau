#pragma once

#include <string>
#include <vector>

namespace jec {

[[noreturn]] void handleError(const std::string& fClass, const std::string& fMessage);
float getFloat(const std::string& token);
unsigned getUnsigned(const std::string& token);
long int getSigned(const std::string& token);
std::string getSection(const std::string& token);
std::vector<std::string> getTokens(const std::string& fLine);
std::string getDefinitions(const std::string& token);

}
