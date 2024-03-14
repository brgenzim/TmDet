#include <string>
#include <vector>
#include <pugixml.hpp>
#include <TmdetXml.hpp>
#include <TmdetStruct.hpp>

using namespace std;
using namespace gemmi;

namespace UniTmp::TmdetLib {

    TmdetXml::TmdetXml() {

    }

    TmdetXml::~TmdetXml() {

    }

    void TmdetXml::read(string path) {

    }

    void TmdetXml::write(string path) {

    }

    TmdetStruct TmdetXml::parse() {
        TmdetStruct tmdet = TmdetStruct((string)"Tm_Alpha",(string)"");

        return tmdet;
    }

    void TmdetXml::compose(TmdetStruct tmdet) {

    }
            
}