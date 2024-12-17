// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <System/FilePaths.hpp>
#include <VOs/SecStrVec.hpp>
#include <VOs/Protein.hpp>

/**
 * @brief namespace for tmdet helpers
 *
 * @namespace Tmdet
 * @namespace Helpers
 */
namespace Tmdet::Helpers {

    /**
     * @brief simple helper class for visualizing 
     *        structure and annotation by pymol
     */
    class Pymol {
        private:
            /**
             * @brief file name of pml file
             */
            std::string pmlFileName;

            /**
             * @brief output stream
             */
            std::ofstream os;

            /**
             * @brief protein value object
             */
            const Tmdet::VOs::Protein& protein;
            
            /**
             * @brief output of the first, header part of the pml file
             * 
             * @param pdbFile 
             */
            void head(std::string pdbFile);

            /**
             * @brief output of the various regions
             */
            void regions();

            /**
             * @brief output secondary structure vectors
             * 
             * @param color 
             */
            void dumpSecStrVec(std::string color);

            /**
             * @brief output of the last part of the pml file
             * 
             */
            void tail();


        public:

            /**
             * @brief Construct a new Pymol object
             * 
             * @param protein 
             */
            explicit Pymol(const Tmdet::VOs::Protein& protein) :
                protein(protein) {}

            /**
             * @brief show the given pml file
             * 
             * @param pdbFile 
             */
            void show(std::string pdbFile);

    };
}