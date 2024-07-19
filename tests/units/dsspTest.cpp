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
#include <Utils/Dssp.hpp>

using namespace std;

string getPath(std::string pdbCode);
void calcDssp(Tmdet::ValueObjects::TmdetStruct& tmdetVO);
bool assertTrue(string testDescription, bool condition, int lineNumber);
vector<string> getResidueNames(gemmi::Chain& chain);

string fileName;

int main() {
    // generic init
    fileName = filesystem::path(__FILE__).filename();
    Tmdet::Services::ConfigurationService::init();

    // Test case 1
    {

        //    │ DSSP compare failed at chain G
        //    │ expected: '---HHH--'
        //    │ actual:   '-HHH----'
        //    │ expected: '---HHH--'
        //    │ actual:   '-HHH----'
        //    │ max expected levenshtein distance (%): 2, actual distance (%): 50
        //    │ max expected levenshtein distance: 0.16, actual distance: 4
        //    │ Failed asserting that false is true.

        // arrange
        string pdbCode = "7ec3";
        string inputPath = getPath(pdbCode);

        gemmi::cif::Document document = gemmi::cif::read(gemmi::MaybeGzipped(inputPath));
        auto pdb = gemmi::make_structure(std::move(document));
        auto tmdetVO = Tmdet::ValueObjects::TmdetStruct(pdb, document);
        tmdetVO.inputPath = inputPath;
        Tmdet::DTOS::TmdetStruct::parse(tmdetVO);

        // action
        calcDssp(tmdetVO);
        // assert
        string expected = "---GGG--";
        auto actual = Tmdet::Utils::Dssp::getDsspOfChain(tmdetVO.chains[3]);
        if (!assertTrue("Verifying 'G' chain of 7ec3", expected == actual, __LINE__)) {
            cout << "           expected: " << expected << endl;
            cout << "        dssp string: " << actual << endl;
        }
    }

    // Test case 2
    {
        // arrange
        string pdbCode = "6e8r";
        string inputPath = getPath(pdbCode);

        gemmi::cif::Document document = gemmi::cif::read(gemmi::MaybeGzipped(inputPath));
        auto pdb = gemmi::make_structure(std::move(document));
        auto tmdetVO = Tmdet::ValueObjects::TmdetStruct(pdb, document);
        tmdetVO.inputPath = inputPath;
        Tmdet::DTOS::TmdetStruct::parse(tmdetVO);

        // action
        calcDssp(tmdetVO);
        // assert
        string expected = "------EEEE---SS---EEEE---";
        auto actual = Tmdet::Utils::Dssp::getDsspOfChain(tmdetVO.chains[3]);
        if (!assertTrue("Verifying 'D' chain of 6e8r", expected == actual, __LINE__)) {
            cout << "           expected: " << expected << endl;
            cout << "        dssp string: " << actual << endl;
        }

        // assert of resiudes (it's an alignment check)
        vector<string> expectedResidues{
            "PRO", "ALA", "ASN", "GLY", "PRO", "ALA", "VAL", "GLN", "PHE", "PHE",
            "LYS", "GLY", "LYS", "ASN", "GLY", "SER", "ALA", "ASP", "GLN", "VAL",
            "ILE", "LEU", "VAL", "THR", "GLN"
        };
        auto actualResidues = getResidueNames(tmdetVO.chains[3].gemmi);
        assertTrue("Verifying residues of 'D' chain of 6e8r", expected == actual, __LINE__);
    }

    // Test case 3
    {
        // arrange
        string pdbCode = "4egy";
        string inputPath = getPath(pdbCode);

        gemmi::cif::Document document = gemmi::cif::read(gemmi::MaybeGzipped(inputPath));
        auto pdb = gemmi::make_structure(std::move(document));
        auto tmdetVO = Tmdet::ValueObjects::TmdetStruct(pdb, document);
        tmdetVO.inputPath = inputPath;
        Tmdet::DTOS::TmdetStruct::parse(tmdetVO);

        // action
        calcDssp(tmdetVO);
        // assert
        string expected = "-------------STGGGGG---HHHHHHHHHHHHHHTTSS-TT-B---HHHHHHHHT--HHHHHHHHHHHHHHTSEEEETTTEEEE-";
        auto actual = Tmdet::Utils::Dssp::getDsspOfChain(tmdetVO.chains[0]);
        if (!assertTrue("Verifying 'A' chain of 4egy", expected == actual, __LINE__)) {
            cout << "           expected: " << expected << endl;
            cout << "        dssp string: " << actual << endl;
        }
    }
    return 0;
}

bool assertTrue(std::string testDescription, bool condition, int lineNumber) {
    std::cout << (condition ? "Passed: " : "Failed: ") << testDescription;
    if (!condition) {
        std::cout << " (at line " << fileName << ":" << lineNumber << ")";
    }
    std::cout << std::endl;

    return condition;
}

void calcDssp(Tmdet::ValueObjects::TmdetStruct& tmdetVO) {
    Tmdet::Utils::Dssp dssp = Tmdet::Utils::Dssp(tmdetVO);
    dssp.calcDsspOnStructure();
    dssp.writeDsspOnStructure();
}

string getPath(std::string pdbCode) {
    Tmdet::Services::ConfigurationService::init();
    auto basePath = Tmdet::Services::ConfigurationService::getValue(Tmdet::Services::ConfigurationService::Keys::PDB_DIRECTORY);
    auto inputPath(basePath);
    inputPath += (string("/") + pdbCode[1] + pdbCode[2]) + "/" + pdbCode + "_updated.cif.gz";

    return inputPath;
}

vector<string> getResidueNames(gemmi::Chain& chain) {
    vector<string> result;

    auto appendAction = [&result](gemmi::Residue& residue) { result.emplace_back(residue.name); };
    for_each(chain.residues.begin(), chain.residues.end(), appendAction);

    return result;
}
