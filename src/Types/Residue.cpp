#include <fstream>
#include <map>
#include <Services/ConfigurationService.hpp>
#include <Services/ChemicalComponentDirectoryService.hpp>
#include <Types/Residue.hpp>

using namespace std;

namespace Tmdet::Types {


    namespace ResidueType {
        map<string, Residue> ChemicalCompoundDictionary;

        Residue getResidue(const string& threeLetterCode) {
            if (Residues.count(threeLetterCode) > 0) {
                return Residues.at(threeLetterCode);
            }
            if (ChemicalCompoundDictionary.count(threeLetterCode) > 0) {
                return ChemicalCompoundDictionary[threeLetterCode];
            }

            Residue residue = Services::ChemicalComponentDirectoryService::getComponentAsResidue(threeLetterCode);
            ChemicalCompoundDictionary[threeLetterCode] = residue;
            return residue;
        }
    };
};

