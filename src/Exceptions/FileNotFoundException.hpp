#ifndef __TMDET_EXCEPTIONS_FILE_NOT_FOUND_EXCEPTION__
#define __TMDET_EXCEPTIONS_FILE_NOT_FOUND_EXCEPTION__

#include <string>
#include <exception>

namespace Tmdet::Exceptions {

    class FileNotFoundException : public std::exception {
        private:
            /**
             * @brief path of the missing file
             * 
             */
            std::string fileName;

        public:
            /**
             * @brief Construct a new File Not Found Exception object
             * 
             * @param fn 
             */
            FileNotFoundException(std::string fn)
                : fileName(fn) {}

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
                std::string message = "File not found: " + fileName;
                return message.c_str();
            }
    };
}

#endif