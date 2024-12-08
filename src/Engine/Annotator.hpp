#pragma once

#include <array>
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <Engine/RegionHandler.hpp>
#include <Types/Region.hpp>
#include <Utils/SecStrVec.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Residue.hpp>
#include <VOs/Membrane.hpp>

/**
 * @brief namespace for tmdet engine
 */
namespace Tmdet::Engine {

    class Annotator {
        private:
            Tmdet::VOs::Protein& protein;
            Tmdet::Engine::RegionHandler regionHandler;

            void run();
            void setChainsType();
            void annotateChains();
            void smoothRegions(Tmdet::VOs::Chain& chain, std::string what);
            void detectLoops(Tmdet::VOs::Chain& chain);
            void detectLoop(Tmdet::VOs::Chain& chain, int beg, int end);
            void detectInterfacialHelices();
            void detectReEntrantLoops(Tmdet::VOs::Chain& chain);
            bool hasHelixLoop(Tmdet::VOs::Chain& chain, int begin, int end);
            void detectTransmembraneHelices(Tmdet::VOs::Chain& chain);
            bool sameSide(Tmdet::VOs::Chain& chain, int beg, int end);
            std::vector<Tmdet::VOs::SecStrVec> getParallelAlphas(Tmdet::VOs::Membrane& membrane);
            bool checkParallel(Tmdet::VOs::SecStrVec& vec, Tmdet::VOs::Membrane& membrane) const;    
            void finalCheck();
            
        public:
            

            explicit Annotator(Tmdet::VOs::Protein& protein) :
                protein(protein),
                regionHandler(Tmdet::Engine::RegionHandler(protein)) {
                    run();
            }
            ~Annotator()=default;
    };
}