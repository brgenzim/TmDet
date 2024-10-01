
#include <string>

#include <config.hpp>
#include <Optim/Organizer.hpp>
#include <Exceptions/NoProteinStructureException.hpp>
#include <Utils/Dssp.hpp>
#include <Utils/Surface.hpp>
#include <Utils/Symmetry.hpp>

namespace Tmdet::Optim {

    void Organizer::main() {
        if (selectChains()) {
            dssp();
            surface();
            checkSymmetry();
            if (!_tmdetVO.tmp) {
                findMembrane();
            }
            if (_tmdetVO.tmp) {
                annotate();
            }
        }
        else {
            throw Tmdet::Exceptions::NoProteinStructureException(_tmdetVO.code);
        }
    }

    unsigned int Organizer::selectChains() {
        unsigned int ret = 0;
        for (auto& chainVO : _tmdetVO.chains) {
            ret += selectChain(chainVO);
        }
        return ret;
    }

    unsigned int Organizer::selectChain(Tmdet::ValueObjects::Chain& chainVO) {
        int nr = 0;
        for (const auto& residueVO : chainVO.residues) {
            nr += (residueVO.hasAllSideChainAtoms()?1:0);
        }
        if (nr < std::stoi(environment.get("TMDET_MIN_NUMBER_OF_RESIDUES_IN_CHAIN",std::to_string(DEFAULT_TMDET_MIN_NUMBER_OF_RESIDUES_IN_CHAIN)))) {
            chainVO.selected = false;
        }
        return chainVO.selected?1:0;
    }

    void Organizer::dssp() {
        Tmdet::Utils::Dssp dssp = Tmdet::Utils::Dssp(_tmdetVO);
    }

    void Organizer::surface() {
        Tmdet::Utils::Surface surf = Tmdet::Utils::Surface(_tmdetVO);
        surf.main();
        surf.setOutsideSurface();
    }

    void Organizer::checkSymmetry() {
        //Tmdet::Utils::Symmetry symmetry;
        //auto result = symmetry.CheckSymmetry(tmdetVO);

    }

    void Organizer::findMembrane() {

    }

    void Organizer::annotate() {
        
    }
}