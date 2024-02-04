#include <unordered_map>
#include <iostream>
#include <string.h>
#include "pdbArgs.hpp"

using namespace std;

namespace UniTmp::PdbLib::Utils {

    /**
     * Constuctor
     */
    pdbArgs::pdbArgs() {
        this->define(false,"h","help","Get list of arguments","bool","false");
    }

    void pdbArgs::define(bool mandatory, string shortFlag, string longFlag,
                string descr, string type, string defaultValue) {
        this->_args.emplace(
            shortFlag.c_str(),
            _pdbArg({mandatory, false, shortFlag, longFlag, descr, 
                type, defaultValue, (defaultValue.length() > 0?defaultValue:"")})
        );
    }

    void pdbArgs::check() {
        bool err = false;
        if (this->_args[(char *)"h"].value == "false") {
            for (const auto& [name, arg] : this->_args){
                if (arg.mandatory && !arg.has) {
                    cerr << "Error: Mandatory argument is missing: " << name << endl;
                    err = true;
                }
            }
        }
        if (err || this->_args[(char *)"h"].value == "true") {
            this->list();
            exit(EXIT_FAILURE);
        }
    }

    void pdbArgs::list() {
        for (const auto& [name, arg] : this->_args){
            cerr << "-" << arg.shortFlag << ", --" << arg.longFlag << (arg.type=="bool"?"":"=") << endl;
            cerr << "\t\t" << arg.description;
            cerr << " [type: '" << arg.type << "']";
            if (arg.mandatory) {
                cerr << " (mandatory)"; 
            }
            if (arg.defaultValue != "") {
                cerr << " (default: " << arg.defaultValue << ")"; 
            }
            cerr << endl;
        }
    }

    void pdbArgs::set(int argc, char *argv[]) {
        char *c=(char *)nullptr;
        char flag[1024];
        string name;
        for(int i=1; i<argc; i++) {
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
                    }
                }
            }
        }
    }

    string pdbArgs::_setValueByLongFlag(char *flag) {
        for(auto [name, arg] : _args) {
            if (arg.longFlag == (string)flag) {
                return name;
            }
        }
        return (char *)nullptr;
    }

    string pdbArgs::_setValueByShortFlag(char *flag) {
        for(auto [name, arg] : _args) {
            if (arg.shortFlag == (string)flag) {
                return name;
            }
        }
        return (char *)nullptr;
    }

    void pdbArgs::_setValue(string name, char *value) {
        if (value != (char *)nullptr) {
            this->_args[name].value = (string)value;
            this->_args[name].has = 1;
        }
        else {
            cerr << "Syntas error in argument list" << endl;
            this->list();
            exit(EXIT_FAILURE);
        }
    }

    bool pdbArgs::getValueAsBool(string name) {
        if (this->_args.contains(name)) {
            if (this->_args[name].type == "bool") {
                return this->_args[name].value == "true"?true:false;
            }
            else {
                cerr << "Argument type error: " << name << endl;
                exit(EXIT_FAILURE);
            }
        }
        cerr << "Argument name error: " << name << endl;
        exit(EXIT_FAILURE);
    }

    int pdbArgs::getValueAsInt(string name) {
        if (this->_args.contains(name)) {
            if (this->_args[name].type == "int") {
                return atoi(this->_args[name].value.c_str());
            }
            else {
                cerr << "Argument type error: " << name << endl;
                exit(EXIT_FAILURE);
            }
        }
        cerr << "Argument name error: " << name << endl;
        exit(EXIT_FAILURE);
    }

    string pdbArgs::getValueAsString(string name) {
        if (this->_args.contains(name)) {
            if (this->_args[name].type == "string") {
                return this->_args[name].value;
            }
            else {
                cerr << "Argument type error: " << name << endl;
                exit(EXIT_FAILURE);
            }
        }
        cerr << "Argument name error: " << name << endl;
        exit(EXIT_FAILURE);
    }
}
