#ifndef __TMDET_UTILS_DSSP__
#define __TMDET_UTILS_DSSP__

#include <array>
#include <string>
#include <any>
#include <gemmi/model.hpp>

#define DSSPTEMP(res) (std::any_cast<dsspTemp&>((res).any.at(dsspTempIdx)))

namespace Tmdet::Utils {

    const std::string dsspTempIdx = "dsspTemp";

    const int t3 = 0;
    const int t4 = 1;
    const int t5 = 2;
    const int ms = 3;
    const int PB = 0;
    const int APB = 1;

    struct dsspTemp {
        gemmi::Residue *residue; //parent residue
        std::array<gemmi::Residue *, 2> to = {nullptr,nullptr}; //first two best hydrogen bond akceptor residue
        std::array<double, 2> energy = {-1e30,-1e30};   //first two best hydrogen bond energy
        std::array<char, 4> ts = {' ',' ',' ',' '};     //t3, t4, t5, ms helpers
        std::array<int, 2> pbs = {-1, -1};              //pb helpers

        explicit dsspTemp(gemmi::Residue *res) :
            residue(res) {};
    };

    void pdbCalcDsspOnStructure(gemmi::Structure& pdb);
    void pdbCalcDsspOnModel(gemmi::Model& model);
    void pdbCalcDsspOnChain(gemmi::Chain& chain);
    void pdbWriteDsspOnStructure(gemmi::Structure& pdb);
    void pdbWriteDsspOnModel(gemmi::Model& model);
    void pdbWriteDsspOnChain(gemmi::Chain& chain);
    void pdbCreateDsspTemp(gemmi::Chain& chain);
    void pdbCreateHydrogenBonds(gemmi::Chain& chain);
    void pdbSetHydrogenBond(gemmi::Residue* donor, gemmi::Residue akceptor, double energy);
    void pdbDetectTurns(gemmi::Chain& chain, int d);
}
#endif