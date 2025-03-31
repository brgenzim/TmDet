// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

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
     * @brief Missing Envirenment Key Exception
     */
    class MissingEnvironmentKeyException : public std::exception {
        private:
            /**
             * @brief name of the missing key
             * 
             */
            std::string key;

        public:

            /**
             * @brief Construct a new Missing Environment Key Exception object
             * 
             * @param k
             */
            MissingEnvironmentKeyException(std::string k)
                : key(k) {}
            
            /**
             * @brief Destroy the Missing Environment Key Exception object
             * 
             */
            ~MissingEnvironmentKeyException()=default;

            
            /**
             * @brief format the error message and throw it
             * 
             * @return const char* 
             */
            const char* what() const throw() {
                std::string message = "Key ("
                        + key + ") not found, nor default value is given";
                return message.c_str();
            }
    };
}
