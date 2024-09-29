#ifndef __TMDET_EXCEPTIONS_MISSING_ENVIRONMENT_KEY_EXCEPTION__
#define __TMDET_EXCEPTIONS_MISSING_ENVIRONMENT_KEY_EXCEPTION__

#include <string>
#include <exception>

namespace Tmdet::Exceptions {

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

#endif