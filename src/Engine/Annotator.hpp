#pragma once

#include <array>
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <Engine/RegionHandler.hpp>
#include <Types/Region.hpp>
#include <Utils/SecStrVec.hpp>
#include <ValueObjects/Protein.hpp>
#include <ValueObjects/Residue.hpp>
#include <ValueObjects/Membrane.hpp>

/**
 * @brief namespace for tmdet engine
 */
namespace Tmdet::Engine {

    class Annotator {
        private:
            Tmdet::ValueObjects::Protein& protein;
            Tmdet::Engine::SideDetector sideDetector;
            Tmdet::Engine::RegionHandler regionHandler;

            void run();
            void setChainsType();
            void annotateChains();
            void smoothRegions(Tmdet::ValueObjects::Chain& chain, std::string what);
            void detectLoops(Tmdet::ValueObjects::Chain& chain);
            void detectLoop(Tmdet::ValueObjects::Chain& chain, int beg, int end);
            void detectInterfacialHelices();
            void detectReEntrantLoops(Tmdet::ValueObjects::Chain& chain);
            bool hasHelixLoop(Tmdet::ValueObjects::Chain& chain, int begin, int end);
            void detectTransmembraneHelices(Tmdet::ValueObjects::Chain& chain);
            //bool isOnOneSide(Tmdet::ValueObjects::Chain& chain, int beg, int end);
            //void extendRegions(Tmdet::Engine::RegionHandler& regionHandler);
            //std::vector<std::array<int,2>> getNumCross(Tmdet::ValueObjects::Chain& chain, int beg, int end);
            bool sameSide(Tmdet::ValueObjects::Chain& chain, int beg, int end);
            //int getTurnResidue(Tmdet::ValueObjects::Chain& chain, int beg, int end);
            //void setTurnResidues(Tmdet::ValueObjects::Chain& chain, std::vector<std::array<int,2>>& numCross);
            //void extendRegion(Tmdet::ValueObjects::Chain& chain, int beg, int end, Tmdet::Engine::RegionHandler& regionHandler);
            std::vector<Tmdet::ValueObjects::SecStrVec> getParallelAlphas(Tmdet::ValueObjects::Membrane& membrane);
            bool checkParallel(Tmdet::ValueObjects::SecStrVec& vec, Tmdet::ValueObjects::Membrane& membrane) const;    
            void finalCheck();
            
        public:
            

            explicit Annotator(Tmdet::ValueObjects::Protein& protein) :
                protein(protein),
                sideDetector(Tmdet::Engine::SideDetector(protein)),
                regionHandler(Tmdet::Engine::RegionHandler(protein)) {
                    run();
            }
            ~Annotator()=default;
    };
}