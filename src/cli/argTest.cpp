#include <iostream>
#include "pdbArgs.hpp"

UniTmp::PdbLib::Utils::pdbArgs setArguments(int argc, char *argv[]);


int main(int argc, char *argv[]) {
    UniTmp::PdbLib::Utils::pdbArgs args = setArguments(argc,argv);
    int i = args.getValueAsInt("a");
    bool b = args.getValueAsBool("b");
    std::cout << "int argument: " << i << std::endl;
    std::cout << "bool argument: " << (b?"true":"false") << std::endl;
}

UniTmp::PdbLib::Utils::pdbArgs setArguments(int argc, char *argv[]) {
    UniTmp::PdbLib::Utils::pdbArgs args;
    args.define(true,"a","a_name","Int test argument","int","");
    args.define(false,"b","b_name","Bool test argument","bool","false");
    args.set(argc,argv);
    args.check();
    return args;
}