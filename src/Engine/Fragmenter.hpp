#pragma once

#include <vector>

#include <gemmi/math.hpp>
#include <System/Arguments.hpp>
#include <VOs/Membrane.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Region.hpp>

/**
 * @brief namespace for tmdet engine
 */
namespace Tmdet::Engine {

    struct _fragmentData {
        unsigned int id;
        unsigned int clusterId;
        bool tmp;
        Tmdet::VOs::Membrane membrane;
        gemmi::Vec3 normal;
        gemmi::Vec3 origo;
        std::vector<Tmdet::VOs::Region> regions;
    };

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

            void runOnFragments(int numFragments);

            void findClusters();

            bool checkAngle(int i, int j);

            int findBestCluster();

            void runOnBestCluster(int bestClusterId);

            void finalize();

            void run();
            void toString();

        public:
            explicit Fragmenter(Tmdet::VOs::Protein& protein, Tmdet::System::Arguments& args) :
                protein(protein),
                args(args) {
                    run();
            }
    };
}
