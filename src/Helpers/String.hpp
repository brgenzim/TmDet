#pragma once

#include <string>
#include <algorithm>
#include <cctype>

namespace Tmdet::Helpers::String {

    extern void toLower(std::string& data);

    extern std::vector<std::string> explode(const std::string& delimiter, std::string source);
}