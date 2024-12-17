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
     * @brief Syntax Error Exception for handing xml files
     */
    class SyntaxErrorException : public std::exception {
        private:
            /**
             * @brief file name where the syntax error present
             * 
             */
            std::string fileName;

            /**
             * @brief line number of the syntax error
             * 
             */
            int lineNum;

            /**
             * @brief description of the syntax error
             * 
             */
            std::string message;
            
        public:

            /**
             * @brief Construct a new Syntax Error Exception object
             * 
             * @param fn 
             * @param ln 
             * @param msg 
             */
            SyntaxErrorException(std::string fn, int ln, std::string msg)
                : fileName(fn),
                  lineNum(ln),
                  message(msg) {}
            
            /**
             * @brief Destroy the Syntax Error Exception object
             * 
             */
            ~SyntaxErrorException()=default;

            
            /**
             * @brief format the error message and throw it
             * 
             * @return const char* 
             */
            const char* what() const throw() {
                std::string fullMessage = (std::string)"Syntax error in file "
                        + fileName
                        + "(line: " + std::to_string(lineNum) + "): "
                        + "Error: " + message;
                return fullMessage.c_str();
            }
    };
}
