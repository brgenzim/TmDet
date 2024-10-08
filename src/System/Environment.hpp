#pragma once

#include <unordered_map>
#include <string>
#include <iostream>

namespace Tmdet::System {
    
    /**
     * @brief handling environment files and variables
     *
     */
    class Environment {
        
        private:
            /**
             * @var key value pairs of environment variables
             * 
             */
            std::unordered_map<std::string, std::string> _envs;

            /**
             * @brief split line into key and variable
             * 
             * @param line 
             * @return std::pair 
             */
            std::pair<std::string,std::string> split(const std::string& line) const;

            /**
             * @brief replace next variable in value
             *        identified by ${...} by the value
             *        of the variable
             * 
             * @param value 
             * @return bool
             */        
            bool nextVariable(std::string& value, const std::string& fileName, const int& lineNum);

            /**
             * @brief read content of the environment 
             *        file and substitutes variables in them
             * 
             * @param envFile 
             */
            void readEnvFile(const std::string& envFile);


        public:
            /**
             * @brief Get the Envs object
             *
             * @throws FileNotFoundException if specified environment file not found
             *         SyntaxErrorException if environment file contains any syntax error
             * 
             * @return std::unordered_map<std::string,std::string> 
             */
            std::unordered_map<std::string,std::string> getEnvs() const {
                return _envs;
            }

            /**
             * @brief initialize environment variables using env file
             *        and shell environment variables
             * 
             * @param envp 
             * @param envFile 
             */
            void init(char** envp, std::string envFile = "");

            /**
             * @brief get value of an environment variable.
             *        If it is not set return default value
             *
             * @throw MissingEnvironmentKeyException if key not found nor default value is given
             * 
             * @param key 
             * @param defaultValue
             * @return std::string 
             */
            std::string get(const std::string& key, std::string defaultValue="");
            
            
            /**
             * @brief set environment variable
             * 
             * @param key 
             * @param value 
             */
            void set(const std::string& key, const std::string& value);

            /**
             * @brief write out all variables and their values
             * 
             * @param os 
             * @param other 
             * @return std::ostream& 
             */
            friend std::ostream& operator<<(std::ostream& os, const Environment& other) {
                for(const auto& [key,value]: other.getEnvs()) {
                    os << key << ": " << value << std::endl;
                }
                return os;
            }

            /**
             * @brief update environment variables
             * 
             * @param envp 
             */
            void updateByEnvVars(char **envp);
    };    
}
