#pragma once

#include <string>
#include <algorithm>
#include <cctype>

namespace Tmdet::Helpers::String {

    static void toLower(std::string& data) {
        std::transform(data.begin(), data.end(), data.begin(),
            [](unsigned char c){ return std::tolower(c); });
    }
}