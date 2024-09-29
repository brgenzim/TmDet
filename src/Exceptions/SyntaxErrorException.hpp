#ifndef __TMDET_EXCEPTIONS_SYNTAX_ERROR_EXCEPTION__
#define __TMDET_EXCEPTIONS_SYNTAX_ERROR_EXCEPTION__

#include <string>
#include <exception>

namespace Tmdet::Exceptions {

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

#endif