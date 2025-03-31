// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <gemmi/model.hpp>
#include <Engine/RegionHandler.hpp>
#include <System/Arguments.hpp>
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
             * @brief protein value object
             */
            Tmdet::VOs::Protein& protein;

            /**
             * @brief command line arguments
             */
            Tmdet::System::Arguments& args;

            /**
             * @brief region handler 
             */
            Tmdet::Engine::RegionHandler& regionHandler;
            
            int numSheets=0;
            std::vector<int> sheetIndex;
            std::vector<std::vector<int>> connectome;
            std::vector<int> numSheetsInBarrels;
            

            /**
             * @brief main entry point of beta annotation
             */
            void run();

            /**
             * @brief collect sheets those are in the membrane
             * 
             */
            void getSheets();

            /**
             * @brief count connected C alpha atoms between sheets
             * 
             */
            void setConnections();

            /**
             * @brief detect connected sheets
             */
            void detectBarrels();

            /**
             * @brief detect sheet that are part of a beta barrel
             * 
             * @param begin 
             * @param sheetNum 
             * @param prevSheet 
             * @param elements 
             * @return int 
             */
            int detectBarrelSheets(int sheetNum, int prevSheet, std::vector<bool>& elements);

            /**
             * @brief Set barrel index to sheets
             * 
             * @param elements 
             */
            void setIndex(std::vector<bool>& elements);

            /**
             * @brief Set elements of a barrel in residue level
             * 
             * @param chain 
             */
            void setBarrel();

            int numConnects(Tmdet::VOs::Chain& chain, int pos);

            /**
             * @brief detect turns between extended regions
             */
            void detectLoops();

            
            
        public:
            /**
             * @brief Construct a new Beta Annotator object
             * 
             * @param chain 
             * @param regionHandler 
             */
            explicit BetaAnnotator(Tmdet::VOs::Protein& protein,
                Tmdet::System::Arguments& args,
                Tmdet::Engine::RegionHandler& regionHandler) :
                protein(protein),
                args(args),
                regionHandler(regionHandler) {
                    run();
                }

            /**
             * @brief Destroy the Beta Annotator object
             */
            ~BetaAnnotator()=default;

            /**
             * @brief detect residues being inside the beta barrel
             */
            void detectBarrelInside(Tmdet::VOs::Chain& chain);
    };
}
