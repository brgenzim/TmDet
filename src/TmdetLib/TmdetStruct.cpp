#include <string>
#include <vector>
#include <iostream>
#include <TmdetStruct.hpp>
#include <TmdetXml.hpp>
#include <TmdetProteinTypes.hpp>

using namespace std;
using namespace gemmi;

namespace UniTmp::TmdetLib {

    TmdetStruct::TmdetStruct() {
    }

    TmdetStruct::~TmdetStruct() {
        
    }

    void TmdetStruct::code(string code) {
        _code = code;
    }
            
    string TmdetStruct::code() {
        return _code;
    }
    
    void TmdetStruct::tmp(bool tmp) {
        _tmp = tmp;
    }

    bool TmdetStruct::tmp() {
        return _tmp;
    }

    void TmdetStruct::read(string path) {
        TmdetXml xml = TmdetXml();
        xml.read(path);
        _tmp = xml.getTmp();
        _code = xml.getCode();
        _date = xml.getCreateDate();
        _modifications = xml.getModifications();
        _qValue = xml.getQvalue();
        std::unordered_map<const char*, TmdetProteinType> tt = TmdetProteinTypes::all;
        for( auto& [key,t] : tt) {
            cerr << ":::" << t.name << "::" << endl;    
        }
        const char* xt = xml.getTmtype().c_str();
        cerr << "+++" << xt << "++" << endl;
        _type = tt[xt];
        cerr << ":::" << tt[xt].name << "::" << endl;
        exit(EXIT_FAILURE);
    }

    void TmdetStruct::write(string path) {
        TmdetXml xml = TmdetXml();
        xml.create();
        xml.setTmp(_tmp);
        xml.setCode(_code);
        xml.setCreateDate(_date);
        xml.setModifications(_modifications);
        xml.setQvalue(_qValue);
        xml.setTmtype(_type.name);
        xml.write(path);
    }
}
