#include <iostream>
#include <string>
#include <filesystem>
#include <memory>
#include <any>
#include <sstream>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>

#include <System/Environment.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <DTOs/TmdetStruct.hpp>
#include <Utils/Dssp.hpp>

using namespace std;

string getPath(std::string pdbCode);
void calcDssp(Tmdet::ValueObjects::TmdetStruct& tmdetVO);
bool assertTrue(string testDescription, bool condition, int lineNumber);
vector<string> getResidueNames(gemmi::Chain& chain);

string fileName;
Tmdet::System::Environment environment;

int main(int argc, char *argv[], char **envp) {
    
    environment.init(envp);
    fileName = filesystem::path(__FILE__).filename();
    
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

        gemmi::Structure pdb; 
        gemmi::cif::Document document;
        auto tmdetVO = Tmdet::ValueObjects::get(inputPath, pdb, document);

        // action
        calcDssp(tmdetVO);
        // assert
        string expected = "---GGG--";
        auto actual = Tmdet::Utils::Dssp::getSecStructAsString(tmdetVO.chains[3]);
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

        gemmi::Structure pdb; 
        gemmi::cif::Document document;
        auto tmdetVO = Tmdet::ValueObjects::get(inputPath, pdb, document);

        // action
        calcDssp(tmdetVO);
        // assert
        string expected = "------EEEE---SS---EEEE---";
        auto actual = Tmdet::Utils::Dssp::getSecStructAsString(tmdetVO.chains[3]);
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

        gemmi::Structure pdb; 
        gemmi::cif::Document document;
        auto tmdetVO = Tmdet::ValueObjects::get(inputPath, pdb, document);

        // action
        calcDssp(tmdetVO);
        // assert
        string expected = "-------------S-TTSTT---HHHHHHHHHHHHHHTTSS-TT-B---HHHHHHHHT--HHHHHHHHHHHHHHTSEEEETTTEEEE-";
        auto actual = Tmdet::Utils::Dssp::getSecStructAsString(tmdetVO.chains[0]);
        if (!assertTrue("Verifying 'A' chain of 4egy", expected == actual, __LINE__)) {
            cout << "           expected: " << expected << endl;
            cout << "        dssp string: " << actual << endl;
        }
    }

    // Test case 4
    {
        // arrange
        string pdbCode = "7e99";
        string inputPath = getPath(pdbCode);

        gemmi::Structure pdb; 
        gemmi::cif::Document document;
        auto tmdetVO = Tmdet::ValueObjects::get(inputPath, pdb, document);

        // action
        calcDssp(tmdetVO);
        // assert
        string expected = "--S-HHHHHHHHHHHHHHS-TTS-HHHHHHHHHHHHHHHHHHH---GGGSGGGG-TT-HHHHHHHHHHHHHHHHHHGGGG-HHHHHHHHHHHHHHHHTSTT--HHHHHHHHHHHHHHHHHHSTT--HHHHHHHHHHHHHHHTTT-";
        auto actual = Tmdet::Utils::Dssp::getSecStructAsString(tmdetVO.chains[3]);
        if (!assertTrue("Verifying 'D' chain of 7e99", expected == actual, __LINE__)) {
            cout << "           expected: " << expected << endl;
            cout << "        dssp string: " << actual << endl;
        }
    }

    // Test case 5
    {
        // arrange
        string pdbCode = "3ee0";
        string inputPath = getPath(pdbCode);

        gemmi::Structure pdb; 
        gemmi::cif::Document document;
        auto tmdetVO = Tmdet::ValueObjects::get(inputPath, pdb, document);

        // action
        calcDssp(tmdetVO);
        // assert
        string expected = "-------STT--TTTGGGT---TTHHHHHHHHH---";
        auto actual = Tmdet::Utils::Dssp::getSecStructAsString(tmdetVO.chains[0]);
        if (!assertTrue("Verifying 'A' chain of 3ee0", expected == actual, __LINE__)) {
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
}

string getPath(std::string pdbCode) {
    
    auto inputPath = environment.get("PDB_DATA_DIR");
    inputPath += (string("/") + pdbCode[1] + pdbCode[2]) + "/" + pdbCode + "_updated.cif.gz";

    return inputPath;
}

vector<string> getResidueNames(gemmi::Chain& chain) {
    vector<string> result;

    auto appendAction = [&result](gemmi::Residue& residue) { result.emplace_back(residue.name); };
    for_each(chain.residues.begin(), chain.residues.end(), appendAction);

    return result;
}

//
// Debug functions
//

string getTempOfChain(Tmdet::ValueObjects::Chain& chain, string key) {
    stringstream result;
    for (Tmdet::ValueObjects::Residue& res : chain.residues) {
        auto& temp = res.temp[key];
        result << ((temp.has_value()) ? any_cast<char>(temp) : ' ');
    }
    return result.str();
}

void printTempsOfChain(Tmdet::ValueObjects::Chain& chain) {
    stringstream result;
    result << "t3: '" << getTempOfChain(chain, "t3") << "'" << endl;
    result << "t4: '" << getTempOfChain(chain, "t4") << "'" << endl;
    result << "t5: '" << getTempOfChain(chain, "t5") << "'" << endl;
    // result << "ss: '" << Tmdet::Utils::Dssp::getSecStructAsString(chain) << "'" << endl;
    cout << result.str();
}

string printHBondsOfChain(Tmdet::ValueObjects::Chain& chain) {
    stringstream result;
    for (Tmdet::ValueObjects::Residue& res : chain.residues) {
        result << res.gemmi.name << "." << res.idx
            << "[" << res.hbond1.toResIdx << ":"  << res.hbond2.toResIdx << "] ";
    }
    result << endl;
    return result.str();
}
