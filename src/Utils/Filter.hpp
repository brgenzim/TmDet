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

    class Filter {
        private:
            Tmdet::VOs::Protein& protein;
            std::string tempDir;
            std::string methodsDir;

        public:
            explicit Filter(Tmdet::VOs::Protein& protein) : 
                protein(protein) {
                    run();
            } ;
            ~Filter()=default;

    };
}
