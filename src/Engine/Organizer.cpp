
#include <string>

#include <Config.hpp>
#include <System/Arguments.hpp>
#include <Engine/Organizer.hpp>
#include <Exceptions/NoProteinStructureException.hpp>
#include <Utils/Dssp.hpp>
#include <Utils/Surface.hpp>
#include <Utils/Symmetry.hpp>
#include <Utils/Oligomer.hpp>

namespace Tmdet::Engine {

    void Organizer::run() {
        logger.debug("Processing Organizer::run()");
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
        logger.debug(" Processed Organizer::run()");
    }

    unsigned int Organizer::selectChains() {
        logger.debug("Processing Organizer::selectChains()");
        unsigned int ret = 0;
        for (auto& chainVO : protein.chains) {
            ret += selectChain(chainVO);
        }
        logger.debug(" Processed Organizer::selectChains(). Return: {}",ret);
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
        logger.debug("Processing Organizer::dssp()");
        auto dssp = Tmdet::Utils::Dssp(protein);
        logger.debug(" Processed Organizer::dssp()");
    }

    void Organizer::surface() {
        logger.debug("Processing Organizer::surface()");
        auto surf = Tmdet::Utils::Surface(protein,args.getValueAsBool("nc"));
        logger.debug(" Processed Organizer::surface()");
    }

    void Organizer::checkSymmetry() {
        logger.debug("Processing Organizer::checkSymmetry()");
        auto oligomerChains = Tmdet::Utils::Oligomer::getHomoOligomerEntities(protein.gemmi);
        if (!oligomerChains.empty()) {
            auto symmetry = Tmdet::Utils::Symmetry(protein);
            
        }
        logger.debug(" Processed Organizer::checkSymmetry()");
    }

    void Organizer::findMembrane() {

    }

    void Organizer::annotate() {
    }
}