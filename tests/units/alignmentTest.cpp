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
            "ILE", "LYS", "ARG", "SER", "LYS", "LYS", "ASN", "SER", "LEU", "ALA",
            "LEU", "SER", "LEU", "THR", "ALA", "ASP", "GLN", "MET", "VAL", "SER",
            "ALA", "LEU", "LEU", "ASP", "ALA", "GLU", "PRO", "PRO", "ILE", "LEU",
            "TYR", "SER", "GLU", "TYR", "ASP", "PRO", "THR", "ARG", "PRO", "PHE",
            "SER", "GLU", "ALA", "SER", "MET", "MET", "GLY", "LEU", "LEU", "THR",
            "ASN", "LEU", "ALA", "ASP", "ARG", "GLU", "LEU", "VAL", "HIS", "MET",
            "ILE", "ASN", "TRP", "ALA", "LYS", "ARG", "VAL", "PRO", "GLY", "PHE",
            "VAL", "ASP", "LEU", "THR", "LEU", "HIS", "ASP", "GLN", "VAL", "HIS",
            "LEU", "LEU", "GLU", "CYS", "ALA", "TRP", "LEU", "GLU", "ILE", "LEU",
            "MET", "ILE", "GLY", "LEU", "VAL", "TRP", "ARG", "SER", "MET", "GLU",
            "HIS", "PRO", "GLY", "LYS", "LEU", "LEU", "PHE", "ALA", "PRO", "ASN",
            "LEU", "LEU", "LEU", "ASP", "ARG", "ASN", "GLN", "GLY", "LYS", "CYS",
            "VAL", "GLU", "GLY", "MET", "VAL", "GLU", "ILE", "PHE", "ASP", "MET",
            "LEU", "LEU", "ALA", "THR", "SER", "SER", "ARG", "PHE", "ARG", "MET",
            "MET", "ASN", "LEU", "GLN", "GLY", "GLU", "GLU", "PHE", "VAL", "CYS",
            "LEU", "LYS", "SER", "ILE", "ILE", "LEU", "LEU", "ASN", "SER", "GLY",
            "VAL", "TYR", "THR", "PHE", "LEU", "SER", "SER", "THR", "LEU", "LYS",
            "SER", "LEU", "GLU", "GLU", "LYS", "ASP", "HIS", "ILE", "HIS", "ARG",
            "VAL", "LEU", "ASP", "LYS", "ILE", "THR", "ASP", "THR", "LEU", "ILE",
            "HIS", "LEU", "MET", "ALA", "LYS", "ALA", "GLY", "LEU", "THR", "LEU",
            "GLN", "GLN", "GLN", "HIS", "GLN", "ARG", "LEU", "ALA", "GLN", "LEU",
            "LEU", "LEU", "ILE", "LEU", "SER", "HIS", "ILE", "ARG", "HIS", "MET",
            "SER", "ASN", "LYS", "GLY", "MET", "GLU", "HIS", "LEU", "TYR", "SER",
            "MET", "LYS", "CYS", "LYS", "ASN", "VAL", "VAL", "PRO", "LEU", "SER",
            "ASP", "LEU", "LEU", "LEU", "GLU", "MET", "LEU", "ASP", "ALA", "HIS",
            "ARG", "LEU", "HIS", "ALA", "PRO", "THR", "SER",
        };
        actual = getResidueNames(tmdetVO.chains[0].gemmi);
        assertTrue("Verifying 'A' chain of 5eit", expected == actual, __LINE__);
    }

    // Test case 2
    {
        auto tmdetVO = createTmdetStruct("6f9w"); // B chain is short enough
        vector<string> expected = {
            "GLY", "PRO", "LEU", "GLY", "SER", "GLY", "LEU", "ALA", "LYS", "TRP",
            "PHE", "GLY", "SER", "ASP", "MSE", "LEU", "GLN", "GLN", "PRO", "LEU",
            "PRO", "SER", "MSE", "PRO", "ALA", "LYS", "VAL", "ILE", "SER", "VAL",
            "ASP", "GLU", "LEU", "GLU", "TYR", "ARG", "GLN",
        };
        auto actual = getResidueNames(tmdetVO.chains[1].gemmi);
        assertTrue("Verifying 'B' chain of 6f9w", expected == actual, __LINE__);
    }

    // Test 3
    {
        auto tmdetVO = createTmdetStruct("7f7g"); // D chain is short
        vector<string> expected = {
            "ACE", "ARG", "ILE", "ARG", "ARG", "ASP", "GLU", "TYR", "LEU", "LYZ",
            "ALA", "ILE", "GLN", "NH2",
        };
        auto actual = getResidueNames(tmdetVO.chains[3].gemmi);
        assertTrue("Verifying 'D' chain of 7f7g", expected == actual, __LINE__);
    }

    // Test 4
    {
        auto tmdetVO = createTmdetStruct("7ec3"); // G chain
        vector<string> expected = {
            "SER", "ASP", "SER", "ASP", "SER", "ASP", "SER", "ASP",
        };
        auto& chain = tmdetVO.chains[3].gemmi;
        auto actual = getResidueNames(chain);
        assertTrue("Verifying 'G' chain of 7ec3", expected == actual, __LINE__);
        // these residues do not have atom lines
        auto predicate = chain.residues[0].atoms.size() == 0
            && chain.residues[1].atoms.size() == 0;
        assertTrue("Verifying last two residues of 'G' chain have no atoms", predicate, __LINE__);
    }

    // Test 5
    {
        auto tmdetVO = createTmdetStruct("4egy"); // A chain
        vector<string> expected = {
            "MET", "HIS", "HIS", "HIS", "HIS", "HIS", "HIS", "LEU", "GLU", "VAL",
            "LEU", "PHE", "GLN", "GLY", "PRO", "LEU", "GLY", "SER", "GLU", "PHE",
            "MET", "LEU", "PRO", "LYS", "TYR", "ALA", "GLN", "VAL", "LYS", "GLU",
            "GLU", "ILE", "SER", "SER", "TRP", "ILE", "ASN", "GLN", "GLY", "LYS",
            "ILE", "LEU", "PRO", "ASP", "GLN", "LYS", "ILE", "PRO", "THR", "GLU",
            "ASN", "GLU", "LEU", "MET", "GLN", "GLN", "PHE", "GLY", "VAL", "SER",
            "ARG", "HIS", "THR", "ILE", "ARG", "LYS", "ALA", "ILE", "GLY", "ASP",
            "LEU", "VAL", "SER", "GLN", "GLY", "LEU", "LEU", "TYR", "SER", "VAL",
            "GLN", "GLY", "GLY", "GLY", "THR", "PHE", "VAL", "ALA",
        };
        auto& chain = tmdetVO.chains[0].gemmi;
        auto actual = getResidueNames(chain);
        assertTrue("Verifying 'A' chain of 4egy", expected == actual, __LINE__);
    }

    // Test 6
    {
        auto tmdetVO = createTmdetStruct("4em2"); // A chain
        vector<string> expected = {
            "THR", "ALA", "ALA", "ALA", "LYS", "PHE", "GLU", "ARG", "GLN", "HIS",
            "MET", "ASP", "SER", "PRO", "ASP", "LEU", "GLY", "THR", "ASP", "ASP",
            "ASP", "ASP", "LYS", "ALA", "MET", "ALA", "ASP", "ILE", "GLY", "SER",
            "ASP", "PHE", "MET", "LEU", "SER", "GLN", "GLU", "PHE", "PHE", "ASN",
            "SER", "PHE", "ILE", "THR", "ILE", "TYR", "ARG", "PRO", "TYR", "LEU",
            "LYS", "LEU", "THR", "GLU", "PRO", "ILE", "LEU", "GLU", "LYS", "HIS",
            "ASN", "ILE", "TYR", "TYR", "GLY", "GLN", "TRP", "LEU", "ILE", "LEU",
            "ARG", "ASP", "ILE", "ALA", "LYS", "HIS", "GLN", "PRO", "THR", "THR",
            "LEU", "ILE", "GLU", "ILE", "SER", "HIS", "ARG", "ARG", "ALA", "ILE",
            "GLU", "LYS", "PRO", "THR", "ALA", "ARG", "LYS", "THR", "LEU", "LYS",
            "ALA", "LEU", "ILE", "GLU", "ASN", "ASP", "LEU", "ILE", "THR", "VAL",
            "GLU", "ASN", "SER", "LEU", "GLU", "ASP", "LYS", "ARG", "GLN", "LYS",
            "PHE", "LEU", "THR", "LEU", "THR", "PRO", "LYS", "GLY", "HIS", "GLU",
            "LEU", "TYR", "GLU", "ILE", "VAL", "CYS", "LEU", "ASP", "VAL", "GLN",
            "LYS", "LEU", "GLN", "GLN", "ALA", "VAL", "VAL", "ALA", "LYS", "THR",
            "ASN", "ILE", "SER", "GLN", "ASP", "GLN", "MET", "GLN", "GLU", "THR",
            "ILE", "ASN", "VAL", "MET", "ASN", "GLN", "ILE", "HIS", "GLU", "ILE",
            "LEU", "LEU", "LYS", "GLU", "ALA", "HIS", "ASN", "ASP",
        };
        auto& chain = tmdetVO.chains[0].gemmi;
        auto actual = getResidueNames(chain);
        assertTrue("Verifying 'A' chain of 4em2", expected == actual, __LINE__);
    }

    // Test 7
    {
        auto tmdetVO = createTmdetStruct("3ee0"); // A chain
        vector<string> expected = {
            "THR", "PHE", "GLY", "SER", "GLY", "GLU", "ALA", "ASP", "CYS", "GLY",
            "LEU", "ARG", "PRO", "LEU", "PHE", "GLU", "LYS", "LYS", "SER", "LEU",
            "GLU", "ASP", "LYS", "THR", "GLU", "ARG", "GLU", "LEU", "LEU", "GLU",
            "SER", "TYR", "ILE", "ASP", "GLY", "ARG",
        };
        auto& chain = tmdetVO.chains[0].gemmi;
        auto actual = getResidueNames(chain);
        assertTrue("Verifying 'A' chain of 3ee0", expected == actual, __LINE__);
    }

    // Test 8
    {
        auto tmdetVO = createTmdetStruct("3f5b"); // A chain
        vector<string> expected = {
            "SER", "ASN", "ALA", "MSE", "MSE", "ILE", "LYS", "ALA", "SER", "THR",
            "ASN", "GLU", "PHE", "ARG", "PHE", "CYS", "PHE", "LYS", "GLN", "MSE",
            "ASN", "LYS", "SER", "GLN", "HIS", "GLU", "LEU", "VAL", "LEU", "GLY",
            "TRP", "ILE", "HIS", "GLN", "PRO", "HIS", "ILE", "ASN", "GLU", "TRP",
            "LEU", "HIS", "GLY", "ASP", "GLY", "LEU", "SER", "ASN", "THR", "ILE",
            "LYS", "ASP", "LEU", "HIS", "GLU", "PHE", "LEU", "ASN", "ASP", "GLY",
            "LYS", "PRO", "TRP", "ALA", "THR", "HIS", "TRP", "ILE", "ALA", "TYR",
            "ASP", "ASN", "GLU", "ILE", "PRO", "PHE", "ALA", "TYR", "LEU", "ILE",
            "THR", "SER", "GLU", "ILE", "GLU", "LYS", "SER", "GLU", "GLU", "TYR",
            "PRO", "ASP", "GLY", "ALA", "VAL", "THR", "LEU", "ASP", "LEU", "PHE",
            "ILE", "CYS", "ARG", "LEU", "ASP", "TYR", "ILE", "GLY", "LYS", "GLY",
            "LEU", "SER", "VAL", "GLN", "MSE", "ILE", "HIS", "GLU", "PHE", "ILE",
            "LEU", "SER", "GLN", "PHE", "SER", "ASP", "THR", "LYS", "ILE", "VAL",
            "LEU", "ILE", "ASN", "PRO", "GLU", "ILE", "SER", "ASN", "GLU", "ARG",
            "ALA", "VAL", "HIS", "VAL", "TYR", "LYS", "LYS", "ALA", "GLY", "PHE",
            "GLU", "ILE", "ILE", "GLY", "GLU", "PHE", "ILE", "ALA", "SER", "TRP",
            "HIS", "PRO", "VAL", "PRO", "HIS", "TYR", "LYS", "MSE", "LYS", "LEU",
            "CYS", "ILE", "GLU", "ASP", "LEU", "LYS", "LYS", "GLN", "ARG", "LEU",
            "SER", "ALA",
        };
        auto& chain = tmdetVO.chains[0].gemmi;
        auto actual = getResidueNames(chain);
        assertTrue("Verifying 'A' chain of 3f5b", expected == actual, __LINE__);
    }

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
    Tmdet::Services::ConfigurationService::init();
    auto basePath = Tmdet::Services::ConfigurationService::getValue(Tmdet::Services::ConfigurationService::Keys::PDB_DIRECTORY);
    auto inputPath(basePath);
    inputPath += (string("/") + pdbCode[1] + pdbCode[2]) + "/" + pdbCode + "_updated.cif.gz";

    gemmi::Structure pdb; 
    gemmi::cif::Document document;
    auto tmdetVO = Tmdet::ValueObjects::get(inputPath, pdb, document);

    return tmdetVO;
}
