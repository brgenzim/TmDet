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
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define STEMP(a,b) (any_cast<surfTemp&>(a.temp["surfTemp"]).b)

    struct surfTemp {
        vector<double> d;
        vector<double> d2;
        vector<double> beta;
        vector<double> arc1;
        vector<double> arc2;
        vector<int> sorted;
    };

    class Surface {
        private:
            Tmdet::ValueObjects::TmdetStruct& tmdetVO;
            void createTempData();
            void calcTempDataOfAtom(Tmdet::ValueObjects::Atom& atom);
            void calcSurfaceOfAtom(Tmdet::ValueObjects::Atom& atom);


        public:
            Surface(Tmdet::ValueObjects::TmdetStruct& _tmdetVO) : tmdetVO(_tmdetVO) {} ;
            ~Surface() {};

            void main();
    };
}
#endif