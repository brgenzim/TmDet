
#include <string>

#include <Config.hpp>
#include <DTOs/Protein.hpp>
#include <Engine/Annotator.hpp>
#include <Engine/Organizer.hpp>
#include <Exceptions/NoProteinStructureException.hpp>
#include <System/Arguments.hpp>
#include <Utils/Dssp.hpp>
#include <Utils/Surface.hpp>
#include <Utils/Symmetry.hpp>
#include <Utils/Oligomer.hpp>

namespace Tmdet::Engine {

    void Organizer::run() {
        DEBUG_LOG("Processing Organizer::run()");
        if (selectChains()) {
            auto optimizer = Optimizer(protein);
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
        DEBUG_LOG(" Processed Organizer::run()");
    }

    unsigned int Organizer::selectChains() {
        DEBUG_LOG("Processing Organizer::selectChains()");
        unsigned int ret = 0;
        for (auto& chain : protein.chains) {
            ret += selectChain(chain);
        }
        DEBUG_LOG(" Processed Organizer::selectChains(). Return: {}",ret);
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
        DEBUG_LOG("Processing Organizer::dssp()");
        auto dssp = Tmdet::Utils::Dssp(protein);
        DEBUG_LOG(" Processed Organizer::dssp()");
    }

    void Organizer::surface() {
        DEBUG_LOG("Processing Organizer::surface()");
        auto surf = Tmdet::Utils::Surface(protein,args.getValueAsBool("nc"));
        DEBUG_LOG(" Processed Organizer::surface()");
    }

    void Organizer::checkSymmetry(Tmdet::Engine::Optimizer& optimizer) {
        DEBUG_LOG("Processing Organizer::checkSymmetry()");
        if (auto oligomerChains = Tmdet::Utils::Oligomer::getHomoOligomerEntities(protein.gemmi); !oligomerChains.empty()) {
            auto symmetry = Tmdet::Utils::Symmetry(protein);
            auto axes = symmetry.getMembraneAxes();
            for(auto& normal: axes) {
                optimizer.setNormal(normal);
                optimizer.testMembraneNormal();
                optimizer.setMembranesToProtein();
            }
        }
        DEBUG_LOG(" Processed Organizer::checkSymmetry()");
    }

    void Organizer::annotate() {
        DEBUG_LOG("Processing Organizer::annotate()");
        Tmdet::DTOs::Protein::transform(protein);
        auto annotator = Tmdet::Engine::Annotator(protein);
        annotator.detectSides();
        //TODO
        annotator.getRegions();
        DEBUG_LOG(" Processed Organizer::annotate()");
    }
}