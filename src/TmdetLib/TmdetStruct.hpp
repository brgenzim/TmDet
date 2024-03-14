#ifndef __UNITMP_TMDETLIB_TMDET_STRUCT__
#define __UNITMP_TMDETLIB_TMDET_STRUCT__

#include <string>
#include <vector>
#include <TmdetRegionTypes.hpp>
#include <TmdetProteinTypes.hpp>
#include <TmdetMembrane.hpp>
#include <TmdetMembraneTypes.hpp>

using namespace std;

namespace UniTmp::TmdetLib {

    struct _tmdetRegion {
        int beg;
        string begi;
        int end;
        string endi;
        int rbeg;
        int rend;
        TmdetRegionType type;
    };

    struct _tmdetChain {
        string id;
        bool selected;
        int numtm;
        vector<_tmdetRegion> regions;
    };

    struct _tmdetModification {
        string date;
        string descr;
    };

    class TmdetStruct {
        private:
            string _code;
            bool _tmp;
            string _date;
            vector<_tmdetModification> _modifications;
            TmdetProteinType _type;
            double _qValue;
            TmdetMembrane _membrane;

        
        public:
            TmdetStruct();
            ~TmdetStruct();
            void code(string code);
            string code();
            void tmp(bool tmp);
            bool tmp();
            void read(string path);
            void write(string path);

    };
}

#endif