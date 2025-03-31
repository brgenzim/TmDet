// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include<string>

/**
 * @brief namespace for tmdet helpers
 *
 * @namespace Tmdet
 * @namespace Helpers
 */
namespace Tmdet::Helpers::Gzip {

    /**
     * @brief compress content of source into destination.
     * @attention This function is created for small CIF files intentionally;
     *            source file must be small enough for available memory.
     */
    void compress(const std::string& source, const std::string& destination);

    /**
     * @brief Gzip and write data into destination file.
     */
    void writeFile(const std::string& destination, const std::string& data);

}
