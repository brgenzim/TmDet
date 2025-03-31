// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <vector>

#include <gemmi/math.hpp>
#include <System/Arguments.hpp>
#include <VOs/Membrane.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Region.hpp>

/**
 * @brief namespace for tmdet engine
 *
 * @namespace Tmdet
 * @namespace Engine
 */
namespace Tmdet::Engine {

    /**
     * @brief temporary container for fragments data
     * 
     */
    struct _fragmentData {
        /**
         * @brief identifier 
         */
        unsigned int id;

        /**
         * @brief cluster identifier
         */
        unsigned int clusterId;

        /**
         * @brief tmdet results for the fragment (yes or no)
         */
        bool tmp;

        /**
         * @brief flag if the fragment is in the final tmdet run
         */
        bool final;

        /**
         * @brief tmdet membrane definition for the fragment (if tmp)
         */
        Tmdet::VOs::Membrane membrane;

        /**
         * @brief normal vector of the membrane calculated for the fragment
         */
        gemmi::Vec3 normal;

        /**
         * @brief origo of the membrane calculated for the fragment
         */
        gemmi::Vec3 origo;

        /**
         * @brief annotated regions for the fragment
         */
        std::vector<Tmdet::VOs::Region> regions;

        /**
         * @brief chain idx of regions
         */
        std::vector<int> regionChainIndexes;
    };

    /**
     * @brief the main class for fragment analysis
     */
    class Fragmenter {
        protected:
            
            /**
             * @brief protein value object
             */
            Tmdet::VOs::Protein& protein;

            /**
             * @brief command line arguments
             */
            Tmdet::System::Arguments& args;

            /**
             * @brief fragment data
             */
            std::vector<_fragmentData> data;

            /**
             * @brief depo for saving state of the protein
             */
            std::vector<gemmi::Vec3> depo;

            /**
             * @brief run tmdet algorithm on fragments
             * 
             * @param numFragments 
             */
            void runOnFragments(int numFragments);

            /**
             * @brief find fragments having same membrane normal
             */
            void findClusters();

            /**
             * @brief check angle of the membrane normals of two fragments
             * 
             * @param i 
             * @param j 
             * @return true 
             * @return false 
             */
            bool checkAngle(int i, int j);

            /**
             * @brief find the best cluster of fragments (having the most
             *        membrane regions)
             * 
             * @return int 
             */
            int findBestCluster();

            /**
             * @brief run tmdet algorithm on the best cluster of fragments
             * 
             * @param bestClusterId 
             */
            void runOnBestCluster(int bestClusterId);

            /**
             * @brief finaliye the results (annotate false positive and negative
             *        membrane regions)
             */
            void finalize();

            /**
             * @brief save the original coordinates of the protein
             */
            void saveState();

            /**
             * @brief restore original coordinates of the protein
             * 
             */
            void restoreState();

            /**
             * @brief Get region type of the fragment
             * 
             * @param fr 
             * @return Tmdet::Types::Region 
             */
            Tmdet::Types::Region getRegionType(int fr);
            
            /**
             * @brief main entry point of the fragmenter class
             */
            void run();

        public:

            /**
             * @brief Construct a new Fragmenter object
             * 
             * @param protein 
             * @param args 
             */
            explicit Fragmenter(Tmdet::VOs::Protein& protein, Tmdet::System::Arguments& args) :
                protein(protein),
                args(args) {
                    run();
            }
    };
}
