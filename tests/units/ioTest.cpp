#include <format>
#include <iostream>
#include <filesystem>
#include <gemmi/cifdoc.hpp>
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

    std::ofstream logFile("build/ioTest.log");
    logger.addStream(logFile);
    logger.setLevel(Tmdet::System::level::debug);

    // Init
    auto code{"1afo"};
    auto protein = createTmdetStruct(code);

    // Test case 1 - CIF write test
    {
        protein.tmatrix.rot = {
            1,  0,  0,
            0,  0, -1,
            0,  1,  0
        };
        Tmdet::DTOs::Protein::writeCif(protein, std::format("/tmp/{}_printed.tr.cif.gz", code));
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
