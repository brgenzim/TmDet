#include <string>
#include <DTOs/Dssp.hpp>
#include <ValueObjects/Chain.hpp>


namespace Tmdet::DTOs {

    std::string Dssp::getSecondaryStructure(const Tmdet::ValueObjects::Chain& chain) {
        std::string ret = "";
        for(const auto& residue: chain.residues ) {
            ret += residue.ss.code;
        }
        return ret;
    }
    
}
