#ifndef __TMDET_UTILS_SURFACE__
#define __TMDET_UTILS_SURFACE__

#include <array>
#include <string>
#include <any>
#include <gemmi/model.hpp>
#include <ValueObjects/TmdetStruct.hpp>

using namespace std;

namespace Tmdet::Utils {

#define SURF_PROBSIZE 1.4
#define SURF_ZSLICE 0.05
#define SURF_DIST 0.1
#define VDW(a) (any_cast<double&>(a.temp["vdw"]))

    struct surfNeighbor {
        Tmdet::ValueObjects::Atom atom;
        double d;
        double d2;
        double beta;
    };

    struct surfTemp {
        vector<surfNeighbor> neighbors;
        vector<double> arc1;
        vector<double> arc2;
        vector<int> sorted;
    };

    class Surface {
        private:
            Tmdet::ValueObjects::TmdetStruct& tmdetVO;
            void initTempData();
            void setContacts();
            void setContactsOfAtom(Tmdet::ValueObjects::Atom& a_atom);
            void setNeighbor(Tmdet::ValueObjects::Atom& a_atom, Tmdet::ValueObjects::Atom& b_atom, surfTemp& st);
            void calcSurfaceOfAtom(Tmdet::ValueObjects::Atom& atom,  surfTemp& st);
            bool calcArcsOfAtom(Tmdet::ValueObjects::Atom& a_atom, surfTemp& st, double z);
            double calcSumArcsOfAtom(Tmdet::ValueObjects::Atom& atom, surfTemp& st, bool ss);
            
        public:
            Surface(Tmdet::ValueObjects::TmdetStruct& _tmdetVO) : tmdetVO(_tmdetVO) {} ;
            ~Surface() {};

            void main();
    };
}
#endif