#pragma once
#include <vector>
#include <gemmi/model.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Residue.hpp>

namespace StructVO = Tmdet::VOs;

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
            StructVO::Protein& proteinVO;

            const double CaDistLimit = 8.0;
            
            
            std::vector<_cr> getNeighbors(const StructVO::Residue& residueVO);
            std::vector<_cr> getCAlphaNetwork();
            std::vector<std::vector<int>> createFragments(const unsigned long size);
            void writeBackFragmentInfoToStructure(std::vector<std::vector<int>> clusters, std::vector<_cr> crs);
            void freeTempValues();
            bool sameChain(gemmi::Chain* chain, int chainIdx);


        public:
            explicit Fragment(StructVO::Protein& proteinVO) : proteinVO(proteinVO) {} ;
            ~Fragment()=default;

            void run();
            
    };
}
