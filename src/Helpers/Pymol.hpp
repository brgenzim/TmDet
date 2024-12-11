#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <System/FilePaths.hpp>
#include <VOs/SecStrVec.hpp>
#include <VOs/Protein.hpp>

namespace Tmdet::Helpers {

    class Pymol {
        private:
            std::string pmlFileName;
            std::ofstream os;
            const Tmdet::VOs::Protein& protein;
            const std::vector<std::string> colors = {
                "limon", //MEMB
                "yellow", //HELIX
                "yellow", //BETA
                "salmon", //SIDE1
                "lightblue", //SIDE2
                "orange", //LOOP
                "green", //IFH
                "violetpurple", //MEMBINS
                "teal", //INTERMEMB
                "red", //INSIDE
                "blue", //OUTSIDE
                "pink", //PERIPLASM
                "red", //ERROR - FP
                "blue", //ERROR - FN
                "grey80" //UNK
            };

            void head(std::string pdbFile);
            void regions();
            void dumpSecStrVec(std::string color);
            void tail();


        public:

            explicit Pymol(const Tmdet::VOs::Protein& protein) :
                protein(protein) {}

            void show(std::string pdbFile);

    };
}