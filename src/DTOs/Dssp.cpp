// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <string>
#include <DTOs/Dssp.hpp>
#include <VOs/Chain.hpp>


namespace Tmdet::DTOs {

    std::string Dssp::getSecondaryStructure(const Tmdet::VOs::Chain& chain) {
        std::string ret = "";
        for(const auto& residue: chain.residues ) {
            ret += residue.ss.code;
        }
        return ret;
    }
    
}
