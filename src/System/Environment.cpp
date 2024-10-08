#include <unordered_map>
#include <iostream>
#include <string>
#include <utility>
#include <fstream>
#include <System/Environment.hpp>
#include <Exceptions/SyntaxErrorException.hpp>
#include <Exceptions/FileNotFoundException.hpp>
#include <Exceptions/MissingEnvironmentKeyException.hpp>

namespace Tmdet::System {

    std::pair<std::string,std::string> Environment::split(const std::string& line) const {
        size_t pos;
        if ((pos = line.find("=")) == std::string::npos) {
            return std::pair<std::string,std::string>{"",""};
        }
        std::string key = line.substr(0,pos);
        std::string value = line.substr(pos+1);
        if (value[0] == '"') {
            value = value.substr(1,value.length()-2);
        }
        return std::pair<std::string,std::string>{key,value};
    }

    bool Environment::nextVariable(std::string& value, const std::string& fileName, const int& lineNum) {
        size_t begin;
        if ((begin = value.find("${")) == std::string::npos) {
            return false;
        }
        size_t end;
        if ((end = value.find("}",begin)) == std::string::npos) {
            throw Tmdet::Exceptions::SyntaxErrorException(fileName, lineNum, "missing '}'");
        }
        std::string variable = value.substr(begin+2,end-begin-2);
        if ( ! _envs.contains(variable)) {
            throw Tmdet::Exceptions::SyntaxErrorException(fileName, lineNum, "unknow variable ("+variable+")");
        }
        std::string newValue = value.substr(0,begin)
            + _envs[variable]
            + value.substr(end+1);
        value = newValue;
        return true;
    }

    /**
     * @brief read content of the environment 
     *        file and substitutes variables in them
     * 
     * @param envFile 
     */
    void Environment::readEnvFile(const std::string& envFile) {
        int lineNum = 1;
        std::ifstream file(envFile);
        if (file.is_open()) {
            std::string line;
            while(std::getline(file, line)) {
                std::pair<std::string,std::string> p = split(line);
                if (p.first != "") {
                    std::string key = p.first;
                    std::string value = p.second;
                    while(nextVariable(value, envFile, lineNum)) {}
                    _envs[key] = value;
                }
                lineNum++;
            }
            file.close();
        }
        else {
            throw Tmdet::Exceptions::FileNotFoundException(envFile);
        }
    }

    std::string Environment::get(const std::string& key, std::string defaultValue) {
        if (_envs.contains(key)) {
            return _envs[key];
        }
        if (defaultValue == "") {
            throw Tmdet::Exceptions::MissingEnvironmentKeyException(key);
        }
        return defaultValue;
    }

    void Environment::set(const std::string& key, const std::string& value) {
        if (_envs.contains(key)) {
            _envs[key] = value;
        }
        else {
            _envs.insert({key,value});
        }
    }

    void Environment::updateByEnvVars(char **envp) {
        for (char **env = envp; *env != 0; env++) {
            std::string shellEnv = *env;
            std::pair<std::string,std::string> p = split(shellEnv);
            if (_envs.contains(p.first)) {
                _envs[p.first] = p.second;
            }
            else if ( p.first.substr(0,6) == "TMDET_") {
                _envs.insert({p.first,p.second});
            }
        }
    }

    void Environment::init(char** envp, std::string envFile) {
        //get environment file content and shell environment variables
        try {
                readEnvFile(envFile);
        }
        catch(std::string& msg) {
            std::cerr << "Environment file not found: " << msg << std::endl;
        }
        updateByEnvVars(envp);
    }

}
