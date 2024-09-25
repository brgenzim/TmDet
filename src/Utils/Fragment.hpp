#ifndef __TMDET_UTILS_FRAGMENT__
#define __TMDET_UTILS_FRAGMENT__

#include <gemmi/model.hpp>
#include <ValueObjects/TmdetStruct.hpp>

#define CA_DIST 8.0
#define MAX_DIST 10000000
#define CUT_LIMIT 50

namespace Tmdet::Utils {

    struct _cr {
        int chain_idx;
        int residue_idx;
    };

    class Fragment {
        private:
            Tmdet::ValueObjects::TmdetStruct& tmdetVO;
            
            std::vector<_cr> getNeighbors(Tmdet::ValueObjects::Residue& residueVO);
            std::vector<_cr> getCAlphaNetwork();
            std::vector<std::vector<int>> createFragments(int size);
            void writeBackFragmentInfoToStructure(std::vector<std::vector<int>> clusters, std::vector<_cr> crs);
            void freeTempValues();
            bool sameChain(gemmi::Chain* chain, int chainIdx);


        public:
            Fragment(Tmdet::ValueObjects::TmdetStruct& _tmdetVO) : tmdetVO(_tmdetVO) {} ;
            ~Fragment()=default;

            void run();
            
    };
}
#endif
