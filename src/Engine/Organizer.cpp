
#include <string>

#include <System/Config.hpp>
#include <System/Arguments.hpp>
#include <Engine/Organizer.hpp>
#include <Exceptions/NoProteinStructureException.hpp>
#include <Utils/Dssp.hpp>
#include <Utils/Surface.hpp>
#include <Utils/Symmetry.hpp>

namespace Tmdet::Engine {

    void Organizer::run() {
        if (selectChains()) {
            dssp();
            surface();
            checkSymmetry();
            if (!protein.tmp) {
                findMembrane();
            }
            if (protein.tmp) {
                annotate();
            }
        }
        else {
            throw Tmdet::Exceptions::NoProteinStructureException(protein.code);
        }
    }

    unsigned int Organizer::selectChains() {
        unsigned int ret = 0;
        for (auto& chainVO : protein.chains) {
            ret += selectChain(chainVO);
        }
        return ret;
    }

    unsigned int Organizer::selectChain(Tmdet::ValueObjects::Chain& chainVO) {
        int nr = 0;
        for (const auto& residueVO : chainVO.residues) {
            nr += (residueVO.hasAllSideChainAtoms()?1:0);
        }
        if (nr < std::stoi(environment.get("TMDET_MIN_NUMBER_OF_RESIDUES_IN_CHAIN",DEFAULT_TMDET_MIN_NUMBER_OF_RESIDUES_IN_CHAIN))) {
            chainVO.selected = false;
        }
        return chainVO.selected?1:0;
    }

    void Organizer::dssp() {
        auto dssp = Tmdet::Utils::Dssp(protein);
    }

    void Organizer::surface() {
        auto surf = Tmdet::Utils::Surface(protein,args.getValueAsBool("nc"));
    }

    void Organizer::checkSymmetry() {
        Tmdet::Utils::Symmetry symmetry;
        auto result = symmetry.CheckSymmetry(protein);
    }

    void Organizer::findMembrane() {

    }

    void Organizer::annotate() {
    }
}