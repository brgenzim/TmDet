#include <string>
#include <vector>
#include <format>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <System/Command.hpp>
#include <System/FilePaths.hpp>
#include <ValueObjects/SecStrVec.hpp>
#include <ValueObjects/Protein.hpp>
#include <Helpers/Pymol.hpp>

namespace Tmdet::Helpers {

    void Pymol::head() {
        //os << "from pymol.CGO_membrane import *\n";
        os << "from pymol.cgo import *\n";
        os << "from pymol import cmd\n";
        os << "import math\n";
        os << "import random\n";
        os << "load " << Tmdet::System::FilePaths::pdbOut(protein.code) << std::endl;
        os << "cmd.do(\'hide all\')\n";
        os << "cmd.do(\'run data/cgo_arrow.py\')\n";
    }

    void Pymol::regions() {
        for (const auto& chain: protein.chains) {
            for (const auto& region: chain.regions) {
                os << std::format("cmd.do(\'color {}, (chain {} and resi {}:{})\')\n",
                        colors[region.type.id],chain.id,region.beg_auth_seq_id,region.end_auth_seq_id);
            }
        }
    }

    void Pymol::tail() {
        os << "cmd.do(\'cartoon automatic\')\n";
        os << "cmd.do(\'show cartoon\')\n";

    }

    void Pymol::dumpSecStrVec(std::string color) {
        int counter = 1;
        for (auto& vector : protein.vectors) {
            os << std::format("cgo_arrow [ {}, {}, {}], [ {}, {}, {}], color={}, name={}{}\n",
                    vector.begin.x, vector.begin.y, vector.begin.z,
                    vector.end.x, vector.end.y, vector.end.z,
                    color, vector.type.name, counter );
            counter++;
        }
    }
    
    void Pymol::show() {
        int pid = getpid();
        pmlFileName = std::format("/tmp/tmdet_{}.pml",pid);
        os.open(pmlFileName, std::ios::binary);
        head();
        regions();
        dumpSecStrVec("cyan");
        tail();
        os.close();
        Tmdet::System::Command::run((std::string)"pymol " + pmlFileName);
    }

}
