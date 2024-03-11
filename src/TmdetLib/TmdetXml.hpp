#ifndef __UNITMP_TMDETLIB_TMDET_XML__
#define __UNITMP_TMDETLIB_TMDET_XML__

#include <string>
#include <vector>
#include <pugixml.hpp>
#include <TmdetStruct.hpp>

using namespace std;
using namespace gemmi;

namespace UniTmp::TmdetLib {

    class TmdetXml {
        private:
            pugi::xml_document _doc;

        public:
            TmdetXml();
            ~TmdetXml();

            void read(string path);
            void write(string path);
            TmdetStruct parse();
            void compose(TmdetStruct tmdet);
            
    };
}

#endif