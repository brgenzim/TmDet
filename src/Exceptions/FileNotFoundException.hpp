// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <format>
#include <string>
#include <exception>

/**
 * @brief namespace for tmdet exceptions
 *
 * @namespace Tmdet
 * @namespace Exceptions
 */
namespace Tmdet::Exceptions {

    /**
     * @brief File Not Found Exception
     */
    class FileNotFoundException : public std::exception {
        private:
            /**
             * @brief path of the missing file
             */
            std::string fileName;

            /**
             * @brief to store description
             */
            std::string message;

        public:
            /**
             * @brief Construct a new File Not Found Exception object
             *
             * @param fn
             */
            FileNotFoundException(std::string fn)
                : fileName(fn), message{std::format("File not found: {}", fn)} {
                }

            /**
             * @brief Destroy the File Not Found Exception object
             * 
             */
            ~FileNotFoundException()=default;

            /**
             * @brief format the error message and throw it
             * 
             * @return const char* 
             */
            const char* what() const throw() {
                return message.c_str();
            }
    };
}
