// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <ios>

/**
 * @brief namespace for tmdet exceptions
 *
 * @namespace Tmdet
 * @namespace Exceptions
 */
namespace Tmdet::Exceptions {

    /**
     * @brief simple input output exception
     */
    class IOException : public std::ios_base::failure {
        public:
            /**
             * @brief Construct a new IOException object
             * 
             * @param message 
             */
            IOException(const std::string& message) : std::ios_base::failure(message) {}

            /**
             * @brief Construct a new IOException object
             * 
             * @param message 
             */
            IOException(const char* message) : std::ios_base::failure(message) {}
    };

}
