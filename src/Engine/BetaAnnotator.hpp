#pragma once

#include <gemmi/model.hpp>
#include <Engine/RegionHandler.hpp>
#include <ValueObjects/Chain.hpp>

/**
 * @brief namespace for tmdet engine
 */
namespace Tmdet::Engine {

    class BetaAnnotator {
        private:
            Tmdet::ValueObjects::Chain& chain;
            Tmdet::Engine::RegionHandler& regionHandler;
            std::vector<int> reIndex;

            void run();
            void init();
            void end();
            void detectBarrel();
            double averageOutSurface(int beg, int end);
            double averageBeta(int beg, int end);
            double averageDirection(int beg, int end);
            int setCluster(int pos, int cluster, int count);
            void detectLoops();
            void detectBarrelInside();
        
        public:
            explicit BetaAnnotator(Tmdet::ValueObjects::Chain& chain, Tmdet::Engine::RegionHandler& regionHandler) :
                chain(chain),
                regionHandler(regionHandler) {
                    run();
                }

            ~BetaAnnotator()=default;
    };
}
