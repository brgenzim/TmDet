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
                "red", //SIDE1
                "blue", //SIDE2
                "orange", //LOOP
                "purpleblue", //IFH
                "pink", //MEMBINS
                "teal", //INTERMEMB
                "grey80", //UNK
                "red", //INSIDE
                "blue" //OUTSIDE
            };

            void head();
            void regions();
            void dumpSecStrVec(std::string color);
            void tail();


        public:

            explicit Pymol(const Tmdet::VOs::Protein& protein) :
                protein(protein) {}

            void show();

    };
}