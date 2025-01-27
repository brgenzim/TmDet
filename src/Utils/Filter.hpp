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
            bool filePutContents(std::string filePath, std::string content);
            bool createTempFasta(Tmdet::VOs::Chain& chain );
            int runPhobius(std::string id);
            int runScampi(std::string id);
            int runTMHMM(std::string id);
            int parsePhobius(std::string results);
            int parseScampi(std::string results);
            int parseTMHMM(std::string results);
            void run(bool applyTmFilter);

        public:
            explicit Filter(Tmdet::VOs::Protein& protein, bool applyTmFilter) : 
                protein(protein) {
                    run(applyTmFilter);
            }
            ~Filter()=default;

            
    };
}
