#include <fstream>
#include <map>
#include <functional>
#include <Services/ChemicalComponentDirectoryService.hpp>
#include <Types/Residue.hpp>

namespace Tmdet::Types {

    namespace ResidueType {
        std::map<std::string, Residue, std::less<>> ChemicalCompoundDictionary;

        Residue getResidue(const std::string& threeLetterCode) {
            if (Residues.contains(threeLetterCode)) {
                return Residues.at(threeLetterCode);
            }
            if (ChemicalCompoundDictionary.contains(threeLetterCode)) {
                return ChemicalCompoundDictionary[threeLetterCode];
            }

            Residue residue = Services::ChemicalComponentDirectoryService::getComponentAsResidue(threeLetterCode);
            ChemicalCompoundDictionary[threeLetterCode] = residue;
            return residue;
        }
    };
};

