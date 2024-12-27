// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once
#include <vector>
#include <gemmi/model.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Chain.hpp>
#include <VOs/Residue.hpp>

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
            Tmdet::VOs::Protein& proteinVO;
            int chIdx;
            
            std::vector<_cr> getNeighbors(const Tmdet::VOs::Residue& residueVO);
            std::vector<_cr> getCAlphaNetwork();
            std::vector<std::vector<int>> createFragments(const unsigned long size);
            void writeBackFragmentInfoToStructure(std::vector<std::vector<int>> clusters, std::vector<_cr> crs);
            void freeTempValues();
            bool sameChain(gemmi::Chain* chain, int chainIdx);
            
        public:
            explicit Fragment(Tmdet::VOs::Protein& proteinVO) : 
                proteinVO(proteinVO) {} ;
            ~Fragment()=default;

            int run();            
    };
}
