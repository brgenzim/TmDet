// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <string>
#include <vector>
#include <format>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <System/Command.hpp>
#include <System/FilePaths.hpp>
#include <VOs/SecStrVec.hpp>
#include <VOs/Protein.hpp>
#include <Helpers/Pymol.hpp>

namespace Tmdet::Helpers {

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

    void Pymol::head(std::string pdbFile) {
        os << "from pymol.cgo import *\n";
        os << "from pymol import cmd\n";
        os << "import math\n";
        os << "import random\n";
        os << "load " << pdbFile << std::endl;
        os << "cmd.do(\'hide all\')\n";
        os << "cmd.do(\'run data/cgo_arrow.py\')\n";
    }

    void Pymol::regions() {
        for (const auto& chain: protein.chains) {
            for (const auto& region: chain.regions) {
                os << std::format("cmd.do(\'color {}, (chain {} and resi {}:{})\')\n",
                        colors[region.type.id],chain.id,
                        region.beg.authId,
                        region.end.authId);
            }
        }
    }

    void Pymol::tail() {
        os << "cmd.do(\'cartoon automatic\')\n";
        os << "cmd.do(\'show cartoon\')\n";
        os << "cmd.do(\'show spheres, (segi TM_)\')\n";
        os << "cmd.do(\'zoom " << protein.code << "_updated_tr\')";
    }

    void Pymol::dumpSecStrVec(std::string color) {
        int counter = 1;
        for (auto& vector : protein.secStrVecs) {
            os << std::format("cgo_arrow [ {}, {}, {}], [ {}, {}, {}], color={}, name={}{}\n",
                    vector.begin.x, vector.begin.y, vector.begin.z,
                    vector.end.x, vector.end.y, vector.end.z,
                    color, vector.type.name, counter );
            counter++;
        }
    }
    
    void Pymol::show(std::string pdbFile) {
        int pid = getpid();
        pmlFileName = std::format("/tmp/tmdet_{}.pml",pid);
        os.open(pmlFileName, std::ios::binary);
        head(pdbFile);
        regions();
        dumpSecStrVec("cyan");
        tail();
        os.close();
        Tmdet::System::Command::run((std::string)"pymol " + pmlFileName);
    }

}
