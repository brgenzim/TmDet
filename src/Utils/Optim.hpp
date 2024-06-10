#ifndef __TMDET_UTILS_OPTIM__
#define __TMDET_UTILS_OPTIM__

#include <array>
#include <string>
#include <any>
#include <gemmi/model.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <ValueObjects/Membrane.hpp>

using namespace std;
using namespace gemmi;
using namespace Tmdet::ValueObjects;

namespace Tmdet::Utils {

    class Optim {
        private:
            Membrane membraneVO;
            TmdetStruct &tmdetVO;

            void setDistances();

        public:
            Optim(TmdetStruct& tmdetVO) : tmdetVO(tmdetVO) {}
            ~Optim() {}


    };
}
#endif