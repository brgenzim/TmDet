// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <algorithm>
#include <cctype>

/**
 * @brief namespace for tmdet helpers
 *
 * @namespace Tmdet
 * @namespace Helpers
 */
namespace Tmdet::Helpers::String {

    /**
     * @brief convert characters to lower case in a string
     * 
     * @param data 
     * @return std::string 
     */
    extern std::string toLower(std::string data);

    /**
     * @brief convert characters to upper case in a string
     * 
     * @param data 
     * @return std::string 
     */
    extern std::string toUpper(std::string data);

    /**
     * @brief cut into parts a string using a delimiter
     * 
     * @param delimiter 
     * @param source 
     * @return std::vector<std::string> 
     */
    extern std::vector<std::string> explode(const std::string& delimiter, std::string source);

    /**
     * @brief format amino acid sequence
     * 
     * @param sequence 
     * @param lineLength 
     * @param strLength 
     * @param preprend 
     * @return std::string 
     */
    extern std::string formatSequence(std::string sequence, int lineLength = 50, int strLength = 10, std::string preprend = "\t\t");
}
