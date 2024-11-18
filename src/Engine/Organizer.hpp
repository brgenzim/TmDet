#pragma once

#include <vector>
#include <map>
#include <set>

#include <Engine/Optimizer.hpp>
#include <System/Arguments.hpp>
#include <ValueObjects/Protein.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/Protein.hpp>

/**
 * @brief name space for tmdet engine
 */
namespace Tmdet::Engine {

    /**
     * @brief the main class of the engine organizing the
     * 
     */
    class Organizer {
        private:

            /**
             * @brief structure and tmdet data containing value object
             */
            Tmdet::ValueObjects::Protein& protein;

            /**
             * @brief command line arguments
             */
            Tmdet::System::Arguments& args;

            /**
             * @brief select chains that contains more than 15 residues
             *        and return the number of these chains
             * 
             * @return unsigned int 
             */
            unsigned int selectChains();

            /**
             * @brief get the number of residues having all side chain atoms
             *        and select the chain if it is larger than the defined 
             *        minimum
             * 
             * @param chainVO 
             * @return unsigned int 
             */
            unsigned int selectChain(Tmdet::ValueObjects::Chain& chainVO);

            /**
             * @brief calculate the secondary structure of selected chains
             */
            void dssp();

            /**
             * @brief calculate the solvent accessible surface and outside surface
             *        of the selected chains
             */
            void surface();

            /**
             * @brief check if the protein contains homo oligomer parts or entirely
             *        is that and the rotational axes for that part and check if
             *        the rotational axes can be the membrane normal
             */
            void checkSymmetry(Tmdet::Engine::Optimizer& optimizer);

            /**
             * @brief if no symmetry in the molecule or it can not be the membrane plane
             *        search for the best membrane plane
             */
            void findMembrane();

            /**
             * @brief run the process
             */
            void run();

            void searchForOneTm();
            
        public:
            /**
             * @brief Construct a new Organizer object
             * 
             * @param protein the protein structure
             * @param args    command line arguments
             */
            explicit Organizer(Tmdet::ValueObjects::Protein& protein, Tmdet::System::Arguments& args) :
                protein(protein),
                args(args) {
                    run();
            }

            /**
             * @brief Destroy the Organizer object
             */
            ~Organizer()=default;
            
    };
}
