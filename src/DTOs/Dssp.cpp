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
