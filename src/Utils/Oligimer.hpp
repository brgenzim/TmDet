#ifndef __TMDET_UTILS_OLIGOMER__
#define __TMDET_UTILS_OLIGOMER__

#include <array>
#include <string>
#include <any>
#include <gemmi/model.hpp>
#include <ValueObjects/TmdetStruct.hpp>

using namespace std;

namespace Tmdet::Utils {

    class Oligomer {
        private:
            Tmdet::ValueObjects::TmdetStruct& tmdetVO;

        public:
            Oligomer(Tmdet::ValueObjects::TmdetStruct& _tmdetVO) :
                tmdetVO(_tmdetVO) {}
            ~Oligomer();
    }
}
#endif