// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt


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
#include <Engine/CurvedOptimizer.hpp>
#include <System/Arguments.hpp>
#include <Types/Protein.hpp>
#include <Utils/Dssp.hpp>
#include <Utils/Surface.hpp>
#include <Utils/Symmetry.hpp>
#include <Utils/Oligomer.hpp>

namespace Tmdet::Engine {

    Organizer::~Organizer() {
        DEBUG_LOG("Destroying Organizer");
    }

    void Organizer::run() {
        DEBUG_LOG("Processing Organizer::run()");
        if (protein.numberOfSelectedChains() > 0) {
            surface();
            if (args.getValueAsBool("cm")) {
                DEBUG_LOG("Curved optimization");
                protein.forceSingleMembrane = true;
                optimizer = std::make_unique<CurvedOptimizer>(protein,args);
            }
            else {
                DEBUG_LOG("Plane optimization");
                optimizer = std::make_unique<PlaneOptimizer>(protein,args);
            }
            if (args.getValueAsBool("fr")) {
                protein.forceSingleMembrane = true;
            }
            else if (!args.getValueAsBool("ns")) {
                checkSymmetry();
            }

            if (!protein.tmp) {
                optimizer->searchForMembraneNormal();
                optimizer->setMembranesToProtein();
            }
            if (bool na = args.getValueAsBool("na"); !na && protein.tmp) {
                protein.transform();
                auto annotator = Tmdet::Engine::Annotator(protein, args);
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
            for(auto& normal: axes) {
                optimizer->setNormal(normal);
                optimizer->clear();
                optimizer->testMembraneNormal();
                if (optimizer->getType() == "Curved") {
                    normal *= -1.0;
                    optimizer->setNormal(normal);
                    optimizer->testMembraneNormal();
                }
                optimizer->setMembranesToProtein();
            }
        }
        DEBUG_LOG(" Processed Organizer::checkSymmetry()");
    }

}
