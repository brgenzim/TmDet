#pragma once

#include <string>
#include <exception>

namespace Tmdet::Exceptions {

    class NoProteinStructureException : public std::exception {
        private:
            /**
             * @brief code of the protein
             * 
             */
            std::string _code;
            
        public:

            /**
             * @brief Construct a new No Protein Structure Exception object
             * 
             * @param code
             */
            explicit NoProteinStructureException(std::string code)
                : _code(code) {}
            
            /**
             * @brief Destroy the No Protein Structure Exception object
             * 
             */
            ~NoProteinStructureException()=default;

            
            /**
             * @brief format the error message and throw it
             * 
             * @return const char* 
             */
            const char* what() const throw() {
                std::string message = _code
                        + " has no protein structure";
                return message.c_str();
            }
    };
}
