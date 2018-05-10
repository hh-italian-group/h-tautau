#pragma once

#ifndef CondFormats_JetMETObjects_Utilities_h
#define CondFormats_JetMETObjects_Utilities_h

#include <stdexcept>

#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>
#include <utility>

namespace
{
  void handleError(const std::string& fClass, const std::string& fMessage)
  {
#ifdef STANDALONE
    std::stringstream sserr;
    sserr<<fClass<<" ERROR: "<<fMessage;
    throw std::runtime_error(sserr.str());
#else
    std::cout << fClass << ", " << fMessage << std::endl;
#endif
  }
  //----------------------------------------------------------------------
  inline float getFloat(const std::string& token)
  {
    char* endptr;
    float result = float(strtod (token.c_str(), &endptr));
    if (endptr == token.c_str())
      {
        std::stringstream sserr;
        sserr<<"can't convert token "<<token<<" to float value";
    handleError("getFloat",sserr.str());
      }
    return result;
  }
  //----------------------------------------------------------------------
  inline unsigned getUnsigned(const std::string& token)
  {
    char* endptr;
    unsigned result = unsigned(strtoul (token.c_str(), &endptr, 0));
    if (endptr == token.c_str())
      {
        std::stringstream sserr;
        sserr<<"can't convert token "<<token<<" to unsigned value";
    handleError("getUnsigned",sserr.str());
      }
    return result;
  }
  inline long int getSigned(const std::string& token)
  {
    char* endptr;
    unsigned result = unsigned(strtol (token.c_str(), &endptr, 0));
    if (endptr == token.c_str())
      {
        std::stringstream sserr;
        sserr<<"can't convert token "<<token<<" to signed value";
    handleError("getSigned",sserr.str());
      }
    return result;
  }
  //----------------------------------------------------------------------
  inline std::string getSection(const std::string& token)
  {
    size_t iFirst = token.find ('[');
    size_t iLast = token.find (']');
    if (iFirst != std::string::npos && iLast != std::string::npos && iFirst < iLast)
      return std::string (token, iFirst+1, iLast-iFirst-1);
    return "";
  }
  //----------------------------------------------------------------------
  inline std::vector<std::string> getTokens(const std::string& fLine)
  {
    std::vector<std::string> tokens;
    std::string currentToken;
    for (unsigned ipos = 0; ipos < fLine.length (); ++ipos)
      {
        char c = fLine[ipos];
        if (c == '#') break; // ignore comments
        else if (c == ' ')
          { // flush current token if any
            if (!currentToken.empty())
              {
            tokens.push_back(currentToken);
            currentToken.clear();
              }
          }
        else
          currentToken += c;
      }
    if (!currentToken.empty()) tokens.push_back(currentToken); // flush end
    return tokens;
  }
    
    //----------------------------------------------------------------------
    inline std::string getDefinitions(const std::string& token)
    {
        size_t iFirst = token.find ('{');
        size_t iLast = token.find ('}');
        if (iFirst != std::string::npos && iLast != std::string::npos && iFirst < iLast)
            return std::string (token, iFirst+1, iLast-iFirst-1);
        return "";
    }
}
#endif
