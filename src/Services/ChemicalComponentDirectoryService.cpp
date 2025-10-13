// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <stdexcept>
#include <cstdlib>

#include <gemmi/cifdoc.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/to_cif.hpp>

#include <Config.hpp>
#include <System/Environment.hpp>
#include <System/Command.hpp>
#include <System/ProgressBar.hpp>
#include <Exceptions/MissingEnvironmentKeyException.hpp>
#include <Services/ChemicalComponentDirectoryService.hpp>
#include <Services/CurlWrapperService.hpp>
#include <Types/Atom.hpp>
#include <Types/Residue.hpp>

namespace fs = std::filesystem;

namespace Tmdet::Services {

    bool ChemicalComponentDirectoryService::isBuilt() {
        std::string basePath;
        try {
            basePath = environment.get("TMDET_CC_DIR",DEFAULT_TMDET_CC_DIR);
            return std::filesystem::exists(std::filesystem::path(basePath + "/Z/Z/ZZZ.cif"));
        }
        catch( const Tmdet::Exceptions::MissingEnvironmentKeyException& exception) {
            std::cout << exception.what() << std::endl;
        }
        return false;
    }

    void ChemicalComponentDirectoryService::fetch() {
        auto url = environment.get("TMDET_CC_URL",DEFAULT_TMDET_CC_URL);
        auto dir = environment.get("TMDET_CC_DIR",DEFAULT_TMDET_CC_DIR);
        std::string cmd = (std::string)"mkdir -p " + dir;
        Tmdet::System::Command::run(cmd);
        auto destination = environment.get("TMDET_CC_FILE",DEFAULT_TMDET_CC_FILE);
        std::cout << "Downloading " << url << " ... " << std::flush;
        auto status = CurlWrapperService::download(url, destination);
        if (status != CurlWrapperService::Status::Ok) {
            std::string message("Downloading '");
            message += url + "' to '" + destination
                + "' failed. Please validate the url and the destination directory.";
            throw std::runtime_error(message);
        }
        std::cout << "done." << std::endl;
    }

    Tmdet::Types::Residue ChemicalComponentDirectoryService::getComponentAsResidue(const std::string& threeLetterCode) {
        std::string chemCompDirectory = environment.get("TMDET_CC_DIR",DEFAULT_TMDET_CC_DIR)
                         + "/"
                         + std::string(1, threeLetterCode[0]);
        if (threeLetterCode.size() >= 2) {
            chemCompDirectory += "/" + std::string(1, threeLetterCode[1]);
        }

        gemmi::cif::Document document = getChemicalComponentDocument(threeLetterCode);
        auto& block = document.blocks[0];

        if (!block.has_mmcif_category("_chem_comp_atom") || !block.has_mmcif_category("_chem_comp")) {
            throw std::runtime_error("Expected _chem_comp_atom or _chem_comp category not found");
        }
        auto oneLetterCode = block.find_value("_chem_comp.one_letter_code");

        Tmdet::Types::Residue residue;
        residue.name = threeLetterCode;
        residue.a1 = oneLetterCode->at(0);

        // Use GEMMI's find function to get a table view of the category
        std::vector<std::string> columns = {
            "atom_id", "alt_atom_id", "type_symbol", "pdbx_aromatic_flag"
        };
        auto atomTable = block.find("_chem_comp_atom.", columns);
        for (const auto& row : atomTable) {
            const std::string& type = row.at(2); // type_symbol
            // ignore hydrogen atoms
            if (type == "H") {
                continue;
            }
            std::string atomId = row.at(0); // atom_id
            std::string altAtomId = row.at(1); // alt_atom_id
            if (atomId[0] == '"') {
                // strip off the quote marks
                atomId = std::string(atomId.begin() + 1, atomId.end() - 1);
            }
            Types::Atom atom;
            if (type == "C") {
                const std::string& aromaticFlag = row.at(3); // pdbx_aromatic_flag
                if (aromaticFlag == "Y" || atomId == "C") {
                    atom = Tmdet::Types::Atoms.at("C_CAR");
                } else {
                    atom = Tmdet::Types::Atoms.at("C_ALI");
                }
            } else if (Tmdet::Types::Atoms.contains(type)) {
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
            residue.atoms[altAtomId] = atomData;
        }
        return residue;
    }

    void ChemicalComponentDirectoryService::build() {
        std::string input = environment.get("TMDET_CC_FILE",DEFAULT_TMDET_CC_FILE);
        std::string destDir = environment.get("TMDET_CC_DIR",DEFAULT_TMDET_CC_DIR);
        std::cout << "Preparing to install Chemical Component Directory ... " << std::flush;
        gemmi::cif::Document doc = gemmi::cif::read(gemmi::MaybeGzipped(input));
        std::cout << "done." << std::endl;

        Tmdet::System::ProgressBar pg;
        pg.setTitle("Installing CCD files: ");
        pg.setNumTicks(doc.blocks.size());
        pg.displayPercentage();
        pg.displayTasksDone();

        for (const auto& block : doc.blocks) {
            std::string cifPath = createDir(destDir, block.name) + "/" + block.name + ".cif";
            writeCif(cifPath.c_str(), block);
            pg.tick();
        }
        pg.end();
    }

    std::string ChemicalComponentDirectoryService::createDir(const std::string& destDir, const std::string& cifName) {
        std::string path(destDir);
        path += "/" + std::string(1, cifName[0]);
        if (cifName.size() >= 2) {
            path += "/" + std::string(1, cifName[1]);
        }
        std::string cmd = std::string("mkdir -p ") + path;
        Tmdet::System::Command::run(cmd);
        return path;
    }

    void ChemicalComponentDirectoryService::writeCif(const char *cifPath, const gemmi::cif::Block& block) {
        std::ofstream os;
        gemmi::cif::WriteOptions wo;
        os.open(cifPath);
        write_cif_block_to_stream(os, block, wo);
        os.close();
    }

    gemmi::cif::Document ChemicalComponentDirectoryService::getChemicalComponentDocument(const std::string& threeLetterCode) {

        std::string chemCompDirectory = environment.get("TMDET_CC_DIR",DEFAULT_TMDET_CC_DIR)
                         + "/"
                         + std::string(1, threeLetterCode[0]);
        if (threeLetterCode.size() >= 2) {
            chemCompDirectory += "/" + std::string(1, threeLetterCode[1]);
        }

        return gemmi::cif::read(gemmi::MaybeGzipped(chemCompDirectory + "/" + threeLetterCode + ".cif"));
    }

    std::vector<std::string> ChemicalComponentDirectoryService::getChemicalComponentInfo(const std::string& threeLetterCode, std::vector<std::string> columns) {

        // collected values will be stored in this vector
        std::vector<std::string> chemCompValues;

        gemmi::cif::Document document = getChemicalComponentDocument(threeLetterCode);
        auto& block = document.blocks[0];

        if (!block.has_mmcif_category("_chem_comp")) {
            throw std::runtime_error("Expected _chem_comp category not found");
        }

        // get the first (and probably the only one) row
        auto table = block.find("_chem_comp.", columns);
        auto chemCompRow = table.one();
        // add values of columns
        for (long unsigned int index = 0; index < columns.size(); index++) {
            chemCompValues.push_back(chemCompRow.at(index));
        }

        return chemCompValues;
    }
}
