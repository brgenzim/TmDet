#include <iostream>
#include <string>
#include <filesystem>
#include <stdexcept>
#include <cstdlib>

#include <gemmi/cifdoc.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>

#include <Services/ConfigurationService.hpp>
#include <Services/ChemicalComponentDirectoryService.hpp>
#include <Services/CurlWrapperService.hpp>
#include <Types/Atom.hpp>
#include <Types/Residue.hpp>

namespace fs = std::filesystem;

namespace Tmdet::Services::ChemicalComponentDirectoryService {

    static void download();
    static void split();

    void build() {
        download();
        split();
    }

    bool isBuilt() {
        std::string lastFile(ConfigurationService::ChemicalComponentDirectory + "/Z/Z/ZZZ.cif");

        return fs::exists(fs::path(lastFile));
    }

    void download() {
        auto url = ConfigurationService::ChemicalComponentDirectoryUrl;
        auto destination = ConfigurationService::ChemicalComponentFile;
        std::cout << "Downloading " << url << "..." << std::endl;
        auto status = CurlWrapperService::download(url, destination);
        if (status != CurlWrapperService::Status::Ok) {
            std::string message("Downloading '");
            message += url + "' to '" + destination
                + "' failed. Please validate the url and the destination directory.";
            throw std::runtime_error(message);
        }
        std::cout << "Download completed." << std::endl;
    }

    void split() {
        std::cout << "Splitting " << ConfigurationService::ChemicalComponentFile << "..." << std::endl;

        std::string cmd(ConfigurationService::FragmentCifExec);
        cmd += " -i " + ConfigurationService::ChemicalComponentFile
            + " -d " + ConfigurationService::ChemicalComponentDirectory
            + " -s > /dev/null 2>&1";
        int exitCode = std::system(cmd.c_str());
        if (exitCode != 0) {
            std::string message(ConfigurationService::AppName);
            message += ": command failed: '" + cmd + "'";
            throw std::runtime_error(message);
        }

        std::cout << "Splitting done." << std::endl;
    }

    Tmdet::Types::Residue getComponentAsResidue(const string& threeLetterCode) {
        if (!isBuilt()) {
            build();
        }
        string chemCompDirectory = ConfigurationService::ChemicalComponentDirectory + "/" + std::string(1, threeLetterCode[0]);
        if (threeLetterCode.size() >= 2) {
            chemCompDirectory += "/" + std::string(1, threeLetterCode[1]);
        }

        gemmi::cif::Document document = gemmi::cif::read(gemmi::MaybeGzipped(chemCompDirectory + "/" + threeLetterCode + ".cif"));
        gemmi::cif::Block& block = document.blocks[0];

        if (!block.has_mmcif_category("_chem_comp_atom") || !block.has_mmcif_category("_chem_comp")) {
            throw runtime_error("Expected _chem_comp_atom or _chem_comp category not found");
        }
        auto oneLetterCode = block.find_value("_chem_comp.one_letter_code");

        Tmdet::Types::Residue residue;
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
            std::string atomId = atomLoop.val(row, atomIdCol);
            if (atomId[0] == '"') {
                // strip off the quote marks
                atomId = std::string(atomId.begin() + 1, atomId.end() - 1);
            }
            Types::Atom atom;
            if (type == "C") {
                const std::string& aromaticFlag = atomLoop.val(row, aromaticFlagCol);
                if (aromaticFlag == "Y" || atomId == "C") {
                    atom = atom = Tmdet::Types::Atoms.at("C_CAR");
                } else {
                    atom = Tmdet::Types::Atoms.at("C_ALI");
                }
            } else if (Tmdet::Types::Atoms.count(type) > 0) {
                atom = Tmdet::Types::Atoms.at(type);
            } else {
                atom = Tmdet::Types::AtomType::UNK;
                atom.name = type;
            }
            Types::AtomData atomData;
            atomData.atom = atom;
            // TODO: these values have to be corrected later (here or elsewhere)
            atomData.mean = 0;
            atomData.sds = 0;
            residue.atoms[atomId] = atomData;
        }
        return residue;
    }

}
