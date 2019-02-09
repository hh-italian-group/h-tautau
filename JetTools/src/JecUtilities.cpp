#include "h-tautau/JetTools/include/JecUtilities.h"

#include <sstream>

namespace jec {

[[noreturn]] void handleError(const std::string& fClass, const std::string& fMessage)
{
    std::stringstream sserr;
    sserr<<fClass<<" ERROR: "<<fMessage;
    throw std::runtime_error(sserr.str());
}

float getFloat(const std::string& token)
{
    char* endptr;
    float result = float(strtod (token.c_str(), &endptr));
    if (endptr == token.c_str()) {
        std::stringstream sserr;
        sserr<<"can't convert token "<<token<<" to float value";
        handleError("getFloat",sserr.str());
    }
    return result;
}

unsigned getUnsigned(const std::string& token)
{
    char* endptr;
    unsigned result = unsigned(strtoul (token.c_str(), &endptr, 0));
    if (endptr == token.c_str()) {
        std::stringstream sserr;
        sserr<<"can't convert token "<<token<<" to unsigned value";
        handleError("getUnsigned",sserr.str());
    }
    return result;
}

long int getSigned(const std::string& token)
{
    char* endptr;
    unsigned result = unsigned(strtol (token.c_str(), &endptr, 0));
    if (endptr == token.c_str()) {
        std::stringstream sserr;
        sserr<<"can't convert token "<<token<<" to signed value";
        handleError("getSigned",sserr.str());
    }
    return result;
}

std::string getSection(const std::string& token)
{
    size_t iFirst = token.find ('[');
    size_t iLast = token.find (']');
    if (iFirst != std::string::npos && iLast != std::string::npos && iFirst < iLast)
    return std::string (token, iFirst+1, iLast-iFirst-1);
    return "";
}

std::vector<std::string> getTokens(const std::string& fLine)
{
    std::vector<std::string> tokens;
    std::string currentToken;
    for (unsigned ipos = 0; ipos < fLine.length (); ++ipos) {
        char c = fLine[ipos];
        if (c == '#') break; // ignore comments
        else if (c == ' ') { // flush current token if any
            if (!currentToken.empty()) {
                tokens.push_back(currentToken);
                currentToken.clear();
            }
        } else
            currentToken += c;
    }
    if (!currentToken.empty()) tokens.push_back(currentToken); // flush end
    return tokens;
}

std::string getDefinitions(const std::string& token)
{
    size_t iFirst = token.find ('{');
    size_t iLast = token.find ('}');
    if (iFirst != std::string::npos && iLast != std::string::npos && iFirst < iLast)
        return std::string (token, iFirst+1, iLast-iFirst-1);
    return "";
}

} // namespace jec
