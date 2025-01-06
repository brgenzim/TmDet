// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <gemmi/model.hpp>
#include <VOs/Protein.hpp>

namespace Tmdet::Utils {

    class MyDssp {
        private:
            Tmdet::VOs::Protein& protein;
            void exec();
            void end();
            void setCO();
            void setAngle();

        public:
            explicit MyDssp(Tmdet::VOs::Protein& protein) : 
                protein(protein) {
                    exec();
            } ;
            ~MyDssp()=default;

    };
}
