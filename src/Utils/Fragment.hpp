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


namespace Tmdet::Utils {

    struct _cr {
        int chainIdx;
        int residueIdx;
    };

    class Fragment {
        private:
            Tmdet::VOs::Protein& protein;
            int numFragments;
            int nr;
            std::vector<std::vector<bool>> contactMap;
            std::vector<_cr> cmIndex;
            
            std::vector<_cr> getNeighbors(const Tmdet::VOs::Residue& residue);
            void setContactMap();
            void createFragments();
            void setFragment(Tmdet::VOs::Residue& residue);
            void freeTempValues();
            bool enableMove(Tmdet::VOs::Residue& from, Tmdet::VOs::Residue& to);
            
        public:
            explicit Fragment(Tmdet::VOs::Protein& protein) : 
                protein(protein) {} ;
            ~Fragment()=default;

            int run();            
    };
}
