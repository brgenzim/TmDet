#include <iostream>
#include <string>
#include <filesystem>
#include <memory>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include <Helpers/Gzip.hpp>
#include <System/Environment.hpp>
#include <System/Logger.hpp>
#include <ValueObjects/Protein.hpp>
#include <DTOs/Protein.hpp>

using namespace std;

Tmdet::ValueObjects::Protein createTmdetStruct(string pdbCode);
void assertTrue(string testDescription, bool condition, int lineNumber);

string fileName;
Tmdet::System::Environment environment;
Tmdet::System::Logger logger;

int main(int argc, char *argv[], char **envp) {

    environment.init(envp, ".env");
    fileName = filesystem::path(__FILE__).filename();

    // Init
    auto protein = createTmdetStruct("1a0s");


    // Test case 1 - CIF write test
    {
        Tmdet::DTOs::Protein::writeCif(protein, "/tmp/1a0s_tr.cif");
    }

    // Test case 2 - gzip test
    {
        Tmdet::DTOs::Protein::writeCif(protein, "/tmp/1a0s_tr.cif.gz");
    }

    return 0;
}

void assertTrue(std::string testDescription, bool condition, int lineNumber) {
    std::cout << (condition ? "Passed: " : "Failed: ") << testDescription;
    if (!condition) {
        std::cout << " (at line " << fileName << ":" << lineNumber << ")";
    }
    std::cout << std::endl;
}

Tmdet::ValueObjects::Protein createTmdetStruct(std::string pdbCode) {

    auto inputPath = environment.get("PDB_CIF_DIR");
    inputPath += (string("/") + pdbCode[1] + pdbCode[2]) + "/" + pdbCode + "_updated.cif.gz";
    return Tmdet::DTOs::Protein::get(inputPath);
}
