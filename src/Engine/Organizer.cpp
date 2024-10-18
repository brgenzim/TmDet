
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
            auto optimizer = Optim(protein);
            dssp();
            surface();
            checkSymmetry(optimizer);
            if (!protein.tmp) {
                optimizer.searchForMembraneNormal();
                optimizer.setMembranesToProtein();
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
        for (auto& chain : protein.chains) {
            ret += selectChain(chain);
        }
        logger.debug(" Processed Organizer::selectChains(). Return: {}",ret);
        return ret;
    }

    unsigned int Organizer::selectChain(Tmdet::ValueObjects::Chain& chain) {
        int nr = 0;
        for (const auto& residue : chain.residues) {
            if (chain.selected) {
                nr += (residue.hasAllSideChainAtoms()?1:0);
            }
        }
        if (nr < std::stoi(environment.get("TMDET_MIN_NUMBER_OF_RESIDUES_IN_CHAIN",DEFAULT_TMDET_MIN_NUMBER_OF_RESIDUES_IN_CHAIN))) {
            chain.selected = false;
        }
        return chain.selected?1:0;
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

    void Organizer::checkSymmetry(Tmdet::Engine::Optim& optimizer) {
        logger.debug("Processing Organizer::checkSymmetry()");
        if (auto oligomerChains = Tmdet::Utils::Oligomer::getHomoOligomerEntities(protein.gemmi); !oligomerChains.empty()) {
            auto symmetry = Tmdet::Utils::Symmetry(protein);
            auto axes = symmetry.getMembraneAxes();
            for(auto& normal: axes) {
                optimizer.setNormal(normal);
                optimizer.testMembraneNormal();
                optimizer.setMembranesToProtein();
            }
        }
        logger.debug(" Processed Organizer::checkSymmetry()");
    }

    void Organizer::annotate() {
        logger.debug("Processing Organizer::annotate()");
        //TODO
        logger.debug(" Processed Organizer::annotate()");
    }
}