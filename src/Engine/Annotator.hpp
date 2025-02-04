// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <array>
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <Engine/RegionHandler.hpp>
#include <System/Arguments.hpp>
#include <Types/Region.hpp>
#include <Utils/SecStrVec.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Residue.hpp>
#include <VOs/Membrane.hpp>

/**
 * @brief namespace for tmdet engine
 *
 * @namespace Tmdet
 * @namespace Engine
 */
namespace Tmdet::Engine {

    /**
     * @brief the main class for region annotation
     */
    class Annotator {
        private:

            /**
             * @brief structure and tmdet data containing protein value object
             */
            Tmdet::VOs::Protein& protein;

            /**
             * @brief command line arguments
             */
            Tmdet::System::Arguments& args;

            /**
             * @brief region handler
             */
            Tmdet::Engine::RegionHandler regionHandler;

            double ifhAngleLimit = 15;
            double loopHelixPart = 0.25;

            /**
             * @brief the main entry point for annotation
             */
            void run();

            /**
             * @brief Set chain's types
             */
            void setChainsType();

            /**
             * @brief annotate chain
             */
            void annotateChains();

            /**
             * @brief remedite small errors in raw region data
             * 
             * @param what 
             */
            void smoothRegions(std::string what);

            /**
             * @brief detect small loops between two membrane segments
             */
            void detectLoops();

            /**
             * @brief detect loop between two membrane segment
             * 
             * @param chain 
             * @param beg 
             * @param end 
             */
            void detectLoop(Tmdet::VOs::Chain& chain, int beg, int end);

            double maxDist(Tmdet::VOs::Chain& chain, int pos, int beg, int end, int dir);

            /**
             * @brief detect if a region has element from the other side
             *        of membrane leaflet
             * 
             * @param chain 
             * @param pos 
             * @param beg 
             * @param end 
             * @return true 
             * @return false 
             */
            bool hasOtherSide(Tmdet::VOs::Chain& chain, int pos, int beg, int end);

            /**
             * @brief detect interfacial helices
             */
            void detectInterfacialHelices();

            /**
             * @brief detect reentrant loops
             * 
             * @param chain 
             */
            void detectReEntrantLoops(Tmdet::VOs::Chain& chain);

            /**
             * @brief check if the region contains both helix and loop
             * 
             * @param chain 
             * @param begin 
             * @param end 
             * @return true 
             * @return false 
             */
            bool hasHelixTurnLoop(Tmdet::VOs::Chain& chain, int begin, int end, int& numHelix);

            /**
             * @brief detect alpha helical transmembrane regions
             * 
             * @param chain 
             */
            void detectTransmembraneHelices(Tmdet::VOs::Chain& chain);

            double helixContent(Tmdet::VOs::Chain& chain, int beg, int end);

            /**
             * @brief check if the whole region is in the same side of the membrane
             * 
             * @param chain 
             * @param beg 
             * @param end 
             * @return true 
             * @return false 
             */
            bool sameSide(Tmdet::VOs::Chain& chain, int beg, int end);

            /**
             * @brief Get alpha helical regions parallel to the membrane plane
             * 
             * @param membrane 
             * @return std::vector<Tmdet::VOs::SecStrVec> 
             */
            std::vector<Tmdet::VOs::SecStrVec> getParallelAlphas(Tmdet::VOs::Membrane& membrane);

            /**
             * @brief check if the vector is parallel to the membrane plane
             * 
             * @param vec 
             * @param membrane 
             * @return true 
             * @return false 
             */
            bool checkParallel(Tmdet::VOs::SecStrVec& vec, Tmdet::VOs::Membrane& membrane) const;

            /**
             * @brief calculate average surface of a region
             * 
             * @param chain 
             * @param beg 
             * @param end 
             * @return double 
             */
            double averageSurface(Tmdet::VOs::Chain& chain, int beg, int end);

            /**
             * @brief calculate hydrophobicity momentum of a helix
             * 
             * @param chain 
             * @param beg 
             * @param end 
             * @return double 
             */
            double hydrophocityMomentum(Tmdet::VOs::Chain& chain, int beg, int end);

            /**
             * @brief made a final check on the annotation
             * 
             */
            void finalCheck();

            /**
             * @brief Set the radius of the membrane in the xy plane
             * 
             */
            void setMembraneSize();
            
        public:
            
            /**
             * @brief Construct a new Annotator object
             * 
             * @param protein 
             */
            explicit Annotator(Tmdet::VOs::Protein& protein, Tmdet::System::Arguments& args) :
                protein(protein),
                args(args),
                regionHandler(Tmdet::Engine::RegionHandler(protein)) {
                    run();
            }

            /**
             * @brief Destroy the Annotator object
             */
            ~Annotator()=default;
    };
}