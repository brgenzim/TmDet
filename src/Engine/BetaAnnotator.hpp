// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <gemmi/model.hpp>
#include <Engine/RegionHandler.hpp>
#include <VOs/Chain.hpp>

/**
 * @brief namespace for tmdet engine
 *
 * @namespace Tmdet
 * @namespace Engine
 */
namespace Tmdet::Engine {

    /**
     * @brief class for annotating beta barrel chains
     */
    class BetaAnnotator {
        private:
            /**
             * @brief chain value object
             */
            Tmdet::VOs::Chain& chain;

            /**
             * @brief region handler 
             */
            Tmdet::Engine::RegionHandler& regionHandler;
            
            /**
             * @brief main entry point of beta annotation
             */
            void run();

            /**
             * @brief initialize temporary data
             */
            void init();

            /**
             * @brief end of beta annotation, destroy temporary data
             */
            void end();

            /**
             * @brief detect beta barrel structure using hydrogen bond
             *        network and water accessible outline surface
             */
            void detectBarrel();

            /**
             * @brief calculate the number of hydrogen bridges to other sheet
             * 
             * @param beg 
             * @param end 
             * @return int 
             */
            int otherConnection(int beg, int end);

            /**
             * @brief calculate of the average otside accessible surface of the
             *        extended region
             * 
             * @param beg 
             * @param end 
             * @return double 
             */
            double averageOutSurface(int beg, int end);

            /**
             * @brief calculate the average extended residue content of region
             * 
             * @param beg 
             * @param end 
             * @return double 
             */
            double averageBeta(int beg, int end);

            /**
             * @brief calculate the average direction in z coordinates of the region
             * 
             * @param beg 
             * @param end 
             * @return double 
             */
            double averageDirection(int beg, int end);

            /**
             * @brief search for connected hydrogen bridge network
             * 
             * @param pos 
             * @param cluster 
             * @param count 
             * @return int 
             */
            int setCluster(int pos, int cluster, int count);

            /**
             * @brief detect turns between extended regions
             */
            void detectLoops();

            /**
             * @brief detect residues being inside the beta barrel
             */
            void detectBarrelInside();
            
        public:
            /**
             * @brief Construct a new Beta Annotator object
             * 
             * @param chain 
             * @param regionHandler 
             */
            explicit BetaAnnotator(Tmdet::VOs::Chain& chain, Tmdet::Engine::RegionHandler& regionHandler) :
                chain(chain),
                regionHandler(regionHandler) {
                    run();
                }

            /**
             * @brief Destroy the Beta Annotator object
             */
            ~BetaAnnotator()=default;
    };
}
