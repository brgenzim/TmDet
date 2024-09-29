#ifndef __TMDET_UTILS_FRAGMENT__
#define __TMDET_UTILS_FRAGMENT__

#include <vector>
#include <gemmi/model.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <ValueObjects/Residue.hpp>

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

            const double CaDistLimit = 8.0;
            
            
            std::vector<_cr> getNeighbors(const Tmdet::ValueObjects::Residue& residueVO);
            std::vector<_cr> getCAlphaNetwork();
            std::vector<std::vector<int>> createFragments(const unsigned long size);
            void writeBackFragmentInfoToStructure(std::vector<std::vector<int>> clusters, std::vector<_cr> crs);
            void freeTempValues();
            bool sameChain(gemmi::Chain* chain, int chainIdx);


        public:
            explicit Fragment(Tmdet::ValueObjects::TmdetStruct& _tmdetVO) : tmdetVO(_tmdetVO) {} ;
            ~Fragment()=default;

            void run();
            
    };
}
#endif
