#include <Types/Atom.hpp>
#include <Types/Residue.hpp>
#include <fstream>
#include <map>
#include <gemmi/cifdoc.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>

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

            const string chemCompDirectory = string("CCD-fragments/") + threeLetterCode[0] + "/" + threeLetterCode[1];

            gemmi::cif::Document document = gemmi::cif::read(gemmi::MaybeGzipped(chemCompDirectory + "/" + threeLetterCode + ".cif"));
            gemmi::cif::Block& block = document.blocks[0];

            if (!block.has_mmcif_category("_chem_comp_atom") || !block.has_mmcif_category("_chem_comp")) {
                throw runtime_error("Expected _chem_comp_atom or _chem_comp category not found");
            }
            auto oneLetterCode = block.find_value("_chem_comp.one_letter_code");

            Residue residue;
            residue.name = threeLetterCode;
            residue.a1 = oneLetterCode->at(0);

            auto atomLoop = block.find_loop_item("_chem_comp_atom.comp_id")->loop;
            int loopLength = atomLoop.length();
            int atomIdCol = atomLoop.find_tag("_chem_comp_atom.atom_id");
            int typeSymbolCol = atomLoop.find_tag("_chem_comp_atom.type_symbol");
            int aromaticFlagCol = atomLoop.find_tag("_chem_comp_atom.pdbx_aromatic_flag");
            for (int row = 0; row < loopLength; row++) {
                const std::string& type = atomLoop.val(row, typeSymbolCol);
                // ignore hydrogen atoms
                if (type == "H") {
                    continue;
                }
                const std::string& atomId = atomLoop.val(row, atomIdCol);
                Atom atom;
                if (type == "C") {
                    const std::string& aromaticFlag = atomLoop.val(row, aromaticFlagCol);
                    if (aromaticFlag == "Y" || atomId == "C") {
                        atom = atom = Tmdet::Types::Atoms.at("C_CAR");
                    } else {
                        atom = Tmdet::Types::Atoms.at("C_ALI");
                    }
                } else {
                    atom = Tmdet::Types::Atoms.at(type);
                }
                AtomData atomData;
                atomData.atom = atom;
                // TODO: these values have to be corrected later (here or elsewhere)
                atomData.mean = 0;
                atomData.sds = 0;
                residue.atoms[atomId] = atomData;
            }

            ChemicalCompoundDictionary[threeLetterCode] = residue;
            return residue;
        }
    };
};

