#ifndef __TMDET_OPTIM_ORGANIZER__
#define __TMDET_OPTIM_ORGANIZER__

#include <vector>
#include <map>
#include <set>

#include <ValueObjects/TmdetStruct.hpp>
#include <ValueObjects/Chain.hpp>

namespace Tmdet::Optim {

    class Organizer {
        private:
            Tmdet::ValueObjects::TmdetStruct& _tmdetVO;

            unsigned int selectChains();
            unsigned int selectChain(Tmdet::ValueObjects::Chain& chainVO);
            void dssp();
            void surface();
            void checkSymmetry();
            void findMembrane();
            void annotate();
            
        public:
            explicit Organizer(Tmdet::ValueObjects::TmdetStruct& tmdetVO) :
                _tmdetVO(tmdetVO) {}
            ~Organizer()=default;

            void main();
            
    };
}

#endif
