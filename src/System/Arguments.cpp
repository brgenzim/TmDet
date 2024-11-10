#include <unordered_map>
#include <iostream>
#include <string.h>
#include <System/Arguments.hpp>

namespace Tmdet::System {

    Arguments::Arguments() {
        this->define(false,"h","help","Get list of arguments","bool","false");
    }

    void Arguments::define(bool mandatory, std::string shortFlag, std::string longFlag,
                std::string descr, std::string type, std::string defaultValue) {
        this->_args.emplace(
            shortFlag.c_str(),
            _arg({mandatory, false, shortFlag, longFlag, descr, 
                type, defaultValue, (defaultValue.length() > 0?defaultValue:"")})
        );
    }

    void Arguments::check() {
        bool err = false;
        if (this->_args[(char *)"h"].value == "false") {
            for (const auto& [name, arg] : this->_args){
                if (arg.mandatory && !arg.has) {
                    std::cerr << "Error: Mandatory argument is missing: " << name << std::endl;
                    err = true;
                }
            }
        }
        if (err || this->_args[(char *)"h"].value == "true") {
            this->list();
            exit(EXIT_FAILURE);
        }
    }

    void Arguments::list() {
        for (const auto& [name, arg] : this->_args){
            std::cerr << "-" << arg.shortFlag << ", --" << arg.longFlag << (arg.type=="bool"?"":"=") << std::endl;
            std::cerr << "\t\t" << arg.description;
            std::cerr << " [type: '" << arg.type << "']";
            if (arg.mandatory) {
                std::cerr << " (mandatory)"; 
            }
            if (arg.defaultValue != "") {
                std::cerr << " (default: " << arg.defaultValue << ")"; 
            }
            std::cerr << std::endl;
        }
    }

    void Arguments::set(int argc, char *argv[]) {
        commandLine = "";
        char *c=(char *)nullptr;
        char flag[1024];
        std::string name;
        for(int i=1; i<argc; i++) {
            commandLine += argv[i];
            commandLine += " ";
            if (!strncmp(argv[i],"--",2)) {
                if ((c=strchr(argv[i],'=')) != (char*)nullptr) {
                    strncpy(flag, &argv[i][2], (int)(c-argv[i])-2);
                    c++;
                }
                else {
                    strcpy(flag, &argv[i][2]);
                }
                if ((name = this->_setValueByLongFlag(flag)) != "") {
                    if (this->_args[name].type == "bool") {
                        this->_setValue(name,(char *)"true");
                    }
                    else {
                        this->_setValue(name,c);
                    }
                }
            }
            else {
                if ((name = this->_setValueByShortFlag(&argv[i][1])) != "") {
                    if (this->_args[name].type == "bool") {
                        this->_setValue(name,(char *)"true");
                    }
                    else {
                        this->_setValue(name,argv[++i]);
                        commandLine += argv[i];
                        commandLine += " ";
                    }
                }
            }
        }
    }

    std::string Arguments::_setValueByLongFlag(char *flag) {
        for(auto [name, arg] : _args) {
            if (arg.longFlag == (std::string)flag) {
                return name;
            }
        }
        return std::string();
    }

    std::string Arguments::_setValueByShortFlag(char *flag) {
        for(auto [name, arg] : _args) {
            if (arg.shortFlag == (std::string)flag) {
                return name;
            }
        }
        return std::string();
    }

    void Arguments::_setValue(std::string name, char *value) {
        if (value != (char *)nullptr) {
            this->_args[name].value = (std::string)value;
            this->_args[name].has = 1;
        }
        else {
            std::cerr << "Syntas error in argument list" << std::endl;
            this->list();
            exit(EXIT_FAILURE);
        }
    }

    bool Arguments::getValueAsBool(std::string name) {
        if (this->_args.contains(name)) {
            if (this->_args[name].type == "bool") {
                return this->_args[name].value == "true"?true:false;
            }
            else {
                std::cerr << "Argument type error: " << name << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        std::cerr << "Argument name error: " << name << std::endl;
        exit(EXIT_FAILURE);
    }

    int Arguments::getValueAsInt(std::string name) {
        if (this->_args.contains(name)) {
            if (this->_args[name].type == "int") {
                return atoi(this->_args[name].value.c_str());
            }
            else {
                std::cerr << "Argument type error: " << name << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        std::cerr << "Argument name error: " << name << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string Arguments::getValueAsString(std::string name) {
        if (this->_args.contains(name)) {
            if (this->_args[name].type == "string") {
                return this->_args[name].value;
            }
            else {
                std::cerr << "Argument type error: " << name << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        std::cerr << "Argument name error: " << name << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string Arguments::getCommandLine() const {
        return commandLine;
    }
}
