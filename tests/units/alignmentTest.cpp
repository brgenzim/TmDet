#include <iostream>
#include <string>
#include <filesystem>
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
            "LYS",
            "HIS",
            "LYS",
            "ILE",
            "LEU",
            "HIS",
            "ARG",
            "LEU",
            "LEU",
            "GLN",
            "ASP",
            "SER",
            "SER",
            "SER",
        };
        auto actual = getResidueNames(tmdetVO.chains[2].gemmi);
        assertTrue("Verifying 'C' chain of 5eit", expected == actual, __LINE__);
    }

    // Test case 2
    // {
    //     auto tmdetVO = createTmdetStruct("6f9w"); // B chain is short enough
    // }

    return 0;
}

vector<string> getResidueNames(gemmi::Chain& chain) {
    vector<string> result;

    auto appendAction = [&result](gemmi::Residue residue) { result.emplace_back(residue.name); };
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
    Tmdet::ValueObjects::TmdetStruct tmdetVO = Tmdet::ValueObjects::TmdetStruct(pdb, document);
    tmdetVO.inputPath = inputPath;
    Tmdet::DTOS::TmdetStruct::parse(tmdetVO);

    return tmdetVO;
}
