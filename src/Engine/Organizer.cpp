
#include <string>
#include <memory>

#include <gemmi/metadata.hpp>

#include <Config.hpp>
#include <DTOs/Dssp.hpp>
#include <DTOs/Protein.hpp>
#include <Engine/Annotator.hpp>
#include <Engine/Organizer.hpp>
#include <Engine/Optimizer.hpp>
#include <Engine/PlaneOptimizer.hpp>
#include <Engine/BlendedOptimizer.hpp>
#include <Exceptions/NoProteinStructureException.hpp>
#include <System/Arguments.hpp>
#include <Types/Protein.hpp>
#include <Utils/Dssp.hpp>
#include <Utils/Surface.hpp>
#include <Utils/Symmetry.hpp>
#include <Utils/Oligomer.hpp>

namespace Tmdet::Engine {

    void Organizer::run() {
        DEBUG_LOG("Processing Organizer::run()");
        if (selectChains()) {
            surface();
            auto ssVec = Tmdet::Utils::SecStrVec(protein);
            
            if (args.getValueAsBool("bm")) {
                DEBUG_LOG("Blended optimization");
                protein.forceSingleMembrane = true;
                optimizer = std::make_unique<BlendedOptimizer>(protein);
            }
            else {
                DEBUG_LOG("Plane optimization");
                optimizer = std::make_unique<PlaneOptimizer>(protein);
            }
            checkSymmetry();

            if (!protein.tmp) {
                searchForOneTm();
            }
            if (!protein.tmp) {
                optimizer->searchForMembraneNormal();
                optimizer->setMembranesToProtein();
            }
            if (bool na = args.getValueAsBool("na"); !na && protein.tmp) {
                protein.transform();
                auto annotator = Tmdet::Engine::Annotator(protein);
            }
        }
        else {
            protein.tmp = false;
            protein.qValue = 0;
            protein.type = Tmdet::Types::ProteinType::NOPROTEIN;
        }
        DEBUG_LOG(" Processed Organizer::run()");
    }

    gemmi::Vec3 Organizer::getBestNormal() const {
        return optimizer->getBestNormal();
    }

    unsigned int Organizer::selectChains() {
        DEBUG_LOG("Processing Organizer::selectChains()");
        unsigned int ret = 0;
        for (auto& chain : protein.chains) {
            if (chain.selected) {
                ret += selectChain(chain);
            }
        }
        DEBUG_LOG(" Processed Organizer::selectChains(). Return: {}",ret);
        return ret;
    }

    unsigned int Organizer::selectChain(Tmdet::VOs::Chain& chain) {
        if (protein.gemmi.entities[chain.entityIdx].polymer_type != gemmi::PolymerType::PeptideL ) {
            chain.selected = false;
            return 0;
        }
        int nr = 0;
        int nb = 0;
        for (const auto& residue : chain.residues) {
            nr += (residue.hasAllSideChainAtoms()?1:0);
            nb += (residue.hasOnlyBackBoneAtoms()?1:0);
        }
        DEBUG_LOG("selectChain: id:{} nr:{} nb:{}",chain.id,nr,nb);
        if (nb > nr) {
            chain.type = Tmdet::Types::ChainType::LOW_RES;
            DEBUG_LOG("Low Resolution chain: {}",chain.id);
        }
        if (nr+nb < std::stoi(environment.get("TMDET_MIN_NUMBER_OF_RESIDUES_IN_CHAIN",DEFAULT_TMDET_MIN_NUMBER_OF_RESIDUES_IN_CHAIN))) {
            chain.selected = false;
        }
        return chain.selected?1:0;
    }

    void Organizer::surface() {
        DEBUG_LOG("Processing Organizer::surface()");
        auto surf = Tmdet::Utils::Surface(protein,args.getValueAsBool("nc"));
        //DEBUG_LOG("{}",Tmdet::DTOs::Protein::toString(protein));
        DEBUG_LOG(" Processed Organizer::surface()");
    }

    void Organizer::checkSymmetry() {
        DEBUG_LOG("Processing Organizer::checkSymmetry()");
        if (auto oligomerChains = Tmdet::Utils::Oligomer::getHomoOligomerEntities(protein.gemmi); !oligomerChains.empty()) {
            auto symmetry = Tmdet::Utils::Symmetry(protein);
            auto axes = symmetry.getMembraneAxes();
            for(const auto& normal: axes) {
                optimizer->setNormal(normal);
                optimizer->clear();
                optimizer->testMembraneNormal();
                optimizer->setMembranesToProtein();
            }
        }
        DEBUG_LOG(" Processed Organizer::checkSymmetry()");
    }

    void Organizer::searchForOneTm() {
        DEBUG_LOG("Processing Organizer::searchForOneTm()");
        //TODO: check if protein contains only onel long alpha helix
        //      then search for membrane normal around 30 degree
        //      of the arrow of the helix
        DEBUG_LOG(" Processed Organizer::searchForOneTm()");
    }

}
