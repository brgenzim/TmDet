#include <iostream>
#include <string>
#include <filesystem>
#include <memory>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include <Services/ConfigurationService.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <DTOs/TmdetStruct.hpp>

using namespace std;

Tmdet::ValueObjects::TmdetStruct createTmdetStruct(string pdbCode);
void assertTrue(string testDescription, bool condition, int lineNumber);
vector<string> getResidueNames(gemmi::Chain& chain);

string fileName;

int main() {
    fileName = filesystem::path(__FILE__).filename();

    // Test case 1
    {
        //auto tmdetVO = createTmdetStruct("1bxw");
        auto tmdetVO = createTmdetStruct("5eit");
        vector<string> expected = {
            "LYS", "HIS", "LYS", "ILE",
            "LEU", "HIS", "ARG", "LEU",
            "LEU", "GLN", "ASP", "SER",
            "SER", "SER",
        };
        auto actual = getResidueNames(tmdetVO.chains[2].gemmi);
        assertTrue("Verifying 'C' chain of 5eit", expected == actual, __LINE__);

        expected = {
            "ILE", "LYS", "ARG", "SER", "LYS",
            "LYS", "ASN", "SER", "LEU", "ALA",
            "LEU", "SER", "LEU", "THR", "ALA",
            "ASP", "GLN", "MET", "VAL", "SER",
            "ALA", "LEU", "LEU", "ASP", "ALA",
            "GLU", "PRO", "PRO", "ILE", "LEU",
            "TYR", "SER", "GLU", "TYR", "ASP",
            "PRO", "THR", "ARG", "PRO", "PHE",
            "SER", "GLU", "ALA", "SER", "MET",
            "MET", "GLY", "LEU", "LEU", "THR",
            "ASN", "LEU", "ALA", "ASP", "ARG",
            "GLU", "LEU", "VAL", "HIS",
            "MET", "ILE", "ASN", "TRP",
            "ALA", "LYS", "ARG", "VAL",
            "PRO", "GLY", "PHE", "VAL",
            "ASP", "LEU", "THR", "LEU",
            "HIS", "ASP", "GLN", "VAL",
            "HIS", "LEU", "LEU", "GLU",
            "CYS", "ALA", "TRP", "LEU",
            "GLU", "ILE", "LEU", "MET",
            "ILE", "GLY", "LEU", "VAL",
            "TRP", "ARG", "SER", "MET",
            "GLU", "HIS", "PRO", "GLY",
            "LYS", "LEU", "LEU", "PHE",
            "ALA", "PRO", "ASN", "LEU",
            "LEU", "LEU", "ASP", "ARG",
            "ASN", "GLN", "GLY", "LYS",
            "CYS", "VAL", "GLU", "GLY",
            "MET", "VAL", "GLU", "ILE",
            "PHE", "ASP", "MET", "LEU",
            "LEU", "ALA", "THR", "SER",
            "SER", "ARG", "PHE", "ARG",
            "MET", "MET", "ASN", "LEU",
            "GLN", "GLY", "GLU", "GLU",
            "PHE", "VAL", "CYS", "LEU",
            "LYS", "SER", "ILE", "ILE",
            "LEU", "LEU", "ASN", "SER",
            "GLY", "VAL", "TYR", "THR",
            "PHE", "LEU", "SER", "SER",
            "THR", "LEU", "LYS", "SER",
            "LEU", "GLU", "GLU", "LYS",
            "ASP", "HIS", "ILE", "HIS",
            "ARG", "VAL", "LEU", "ASP",
            "LYS", "ILE", "THR", "ASP",
            "THR", "LEU", "ILE", "HIS",
            "LEU", "MET", "ALA", "LYS",
            "ALA", "GLY", "LEU", "THR",
            "LEU", "GLN", "GLN", "GLN",
            "HIS", "GLN", "ARG", "LEU",
            "ALA", "GLN", "LEU", "LEU",
            "LEU", "ILE", "LEU", "SER",
            "HIS", "ILE", "ARG", "HIS",
            "MET", "SER", "ASN", "LYS",
            "GLY", "MET", "GLU", "HIS",
            "LEU", "TYR", "SER", "MET",
            "LYS", "CYS", "LYS", "ASN",
            "VAL", "VAL", "PRO", "LEU",
            "SER", "ASP", "LEU", "LEU",
            "LEU", "GLU", "MET", "LEU",
            "ASP", "ALA", "HIS", "ARG",
            "LEU", "HIS", "ALA", "PRO",
            "THR", "SER",
        };
        actual = getResidueNames(tmdetVO.chains[0].gemmi);
        assertTrue("Verifying 'A' chain of 5eit", expected == actual, __LINE__);
    }

    // Test case 2
    // {
    //     auto tmdetVO = createTmdetStruct("6f9w"); // B chain is short enough
    // }

    return 0;
}

vector<string> getResidueNames(gemmi::Chain& chain) {
    vector<string> result;

    auto appendAction = [&result](gemmi::Residue& residue) { result.emplace_back(residue.name); };
    for_each(chain.residues.begin(), chain.residues.end(), appendAction);

    return result;
}

void assertTrue(std::string testDescription, bool condition, int lineNumber) {
    std::cout << (condition ? "Passed: " : "Failed: ") << testDescription;
    if (!condition) {
        std::cout << " (at line " << fileName << ":" << lineNumber << ")";
    }
    std::cout << std::endl;
}

Tmdet::ValueObjects::TmdetStruct createTmdetStruct(std::string pdbCode) {
    gemmi::Structure pdb;
    Tmdet::Services::ConfigurationService::init();
    auto basePath = Tmdet::Services::ConfigurationService::getValue(Tmdet::Services::ConfigurationService::Keys::PDB_DIRECTORY);
    auto inputPath(basePath);
    inputPath += (string("/") + pdbCode[1] + pdbCode[2]) + "/" + pdbCode + "_updated.cif.gz";

    gemmi::cif::Document document = gemmi::cif::read(gemmi::MaybeGzipped(inputPath));
    pdb = gemmi::make_structure(std::move(document));
    auto tmdetVO = Tmdet::ValueObjects::TmdetStruct(pdb, document);
    tmdetVO.inputPath = inputPath;
    Tmdet::DTOS::TmdetStruct::parse(tmdetVO);

    return tmdetVO;
}
