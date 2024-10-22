#include <iostream>
#include <string>
#include <filesystem>
#include <memory>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include <System/Environment.hpp>
#include <System/Logger.hpp>
#include <ValueObjects/Protein.hpp>
#include <DTOs/Protein.hpp>

using namespace std;

Tmdet::ValueObjects::Protein createTmdetStruct(string pdbCode);
void assertTrue(string testDescription, bool condition, int lineNumber);
vector<string> getResidueNames(gemmi::Chain& chain);

string fileName;
Tmdet::System::Environment environment;
Tmdet::System::Logger logger;

int main(int argc, char *argv[], char **envp) {

    logger.setLevel(Tmdet::System::level::warn);
    logger.addStream(std::cout);

    environment.init(envp, ".env");
    fileName = filesystem::path(__FILE__).filename();

    // Test case 1
    {
        auto protein = createTmdetStruct("5uow");
        Tmdet::DTOs::Protein::unselectPolymers(protein);
        assertTrue("Verifying chain D of 5uow" , protein.chains[3].selected, __LINE__);
        assertTrue("Verifying chain F of 5uow" , protein.chains[4].selected == false, __LINE__);
        assertTrue("Verifying chain G of 5uow" , protein.chains[5].selected == false, __LINE__);
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
