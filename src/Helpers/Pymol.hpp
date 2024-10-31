#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <System/FilePaths.hpp>
#include <ValueObjects/SecStrVec.hpp>
#include <ValueObjects/Protein.hpp>

namespace Tmdet::Helpers {

    class Pymol {
        private:
            std::string pmlFileName;
            std::ofstream os;
            const Tmdet::ValueObjects::Protein& protein;
            const std::vector<std::string> colors = {
                "limon", //MEMB
                "yellow", //HELIX
                "yellow", //BETA
                "red", //SIDE1
                "blue", //SIDE2
                "orange", //LOOP
                "green", //IFH
                "pink", //MEMBINS
                "teal", //INTERMEMB
                "grey80", //UNK
            };

            void head();
            void regions();
            void dumpSecStrVec(std::string color);
            void tail();


        public:

            explicit Pymol(const Tmdet::ValueObjects::Protein& protein) :
                protein(protein) {}

            void show();

    };
}