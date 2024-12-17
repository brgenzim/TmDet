// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <vector>
#include <map>
#include <set>
#include <memory>

#include <Engine/Optimizer.hpp>
#include <System/Arguments.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Chain.hpp>
#include <VOs/Protein.hpp>

/**
 * @brief namespace for tmdet engine
 *
 * @namespace Tmdet
 * @namespace Engine
 */
namespace Tmdet::Engine {

    /**
     * @brief the main class that organize the whole annotation
     */
    class Organizer {
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
             * @brief pointer to the optimizer object (curved or plain)
             */
            std::unique_ptr<Tmdet::Engine::Optimizer> optimizer;

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
            unsigned int selectChain(Tmdet::VOs::Chain& chainVO);

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
            void checkSymmetry();

            /**
             * @brief run the process
             */
            void run();
            
        public:
            /**
             * @brief Construct a new Organizer object
             * 
             * @param protein the protein structure
             * @param args    command line arguments
             */
            explicit Organizer(Tmdet::VOs::Protein& protein, Tmdet::System::Arguments& args) :
                protein(protein),
                args(args) {
                    run();
            }

            /**
             * @brief Destroy the Organizer object
             */
            ~Organizer();
            
            /**
             * @brief get the normal vector for the best qValue
             */
            gemmi::Vec3 getBestNormal() const;

    };
}
