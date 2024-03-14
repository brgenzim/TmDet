#include <string>
#include <vector>
#include <TmdetStruct.hpp>

using namespace std;
using namespace gemmi;

namespace UniTmp::TmdetLib {

    TmdetStruct::TmdetStruct(string protType, string membType) {
        std::unordered_map<const char*, TmdetProteinType> types = TmdetProteinTypes::all;
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

}
