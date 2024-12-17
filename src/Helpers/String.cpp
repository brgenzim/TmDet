// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <string>
#include <vector>
#include <algorithm>
#include <cctype>

namespace Tmdet::Helpers::String {

    std::string toLower(std::string data) {
        std::transform(data.begin(), data.end(), data.begin(),
            [](unsigned char c){ return std::tolower(c); });
        return data;
    }

    std::string toUpper(std::string data) {
        std::transform(data.begin(), data.end(), data.begin(),
            [](unsigned char c){ return std::toupper(c); });
        return data;
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

    std::string formatSequence(std::string sequence, int lineLength, int strLength, std::string preprend) {
        std::string ret = preprend;
        for (unsigned long int i = 0; i<sequence.length(); i++) {
            ret += sequence[i];
            if ((i+1)%lineLength == 0) {
                ret += '\n';
                ret += preprend;
            }
            else if ((i+1)%strLength == 0) {
                ret += ' ';
            }
        }
        return ret;
    }
}