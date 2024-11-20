#pragma once

#include <string>
#include <algorithm>
#include <cctype>

namespace Tmdet::Helpers::String {

    extern std::string toLower(std::string data);

    extern std::string toUpper(std::string data);

    extern std::vector<std::string> explode(const std::string& delimiter, std::string source);

    extern std::string formatSequence(std::string sequence, int lineLength = 50, int strLength = 10);
}