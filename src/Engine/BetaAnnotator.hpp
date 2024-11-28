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

            void run();
            void init();
            void end();
            void detectBarrel();
            int setCluster(int pos, int cluster, int count);
            bool sameDirection(int p1, int p2);
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
