#include <string>
#include <vector>
#include <algorithm>
#include <cctype>

namespace Tmdet::Helpers::String {

    void toLower(std::string& data) {
        std::transform(data.begin(), data.end(), data.begin(),
            [](unsigned char c){ return std::tolower(c); });
    }

    std::vector<std::string> explode(const std::string& delimiter, std::string source) {
        std::vector<std::string> tokens;
        size_t pos = 0;
        std::string token;
        while ((pos = source.find(delimiter)) != std::string::npos) {
            token = source.substr(0, pos);
            tokens.push_back(token);
            source.erase(0, pos + delimiter.length());
        }
        tokens.push_back(source);

        return tokens;
    }
}