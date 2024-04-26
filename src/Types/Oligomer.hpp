#ifndef __TMDET_TYPES_OLIGOMER__
#define __TMDET_TYPES_OLIGOMER__

#include <unordered_map>
#include <string>

using namespace std;

namespace Tmdet::Types {

    struct Oligomer {
        string name;
        string description;
    };

    namespace OligomerType {
        const Oligomer MONOMER = {
            "Monomer", 
            "There is only one chain in the protein"
        };
        const Oligomer HOMO_OLIGOMER = {
            "HomoOligomer",
            "There is only one chain type but several times"
        };
        const Oligomer HETERO_OLIGOMER = {
            "HeteroOligomer",
            "There several different chains in the protein"
        };
        const Oligomer HOMO_HETERO_OLIGOMER = {
            "HomoHeteroOligomer",
            "There several different chains in the protein but multiple times"
        };
    }

    const unordered_map<string, const Oligomer> Oligomers = {
        { "Monomer", OligomerType::MONOMER },
        { "HomoOligomer", OligomerType::HOMO_OLIGOMER },
        { "HeteroOligomer", OligomerType::HETERO_OLIGOMER },
        { "HomoHeteroOligomer", OligomerType::HOMO_HETERO_OLIGOMER }
    };

}

#endif
