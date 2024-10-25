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
void printDocument(std::ostream& outputStream, const std::string&);

string fileName;
Tmdet::System::Environment environment;
Tmdet::System::Logger logger;

int main(int argc, char *argv[], char **envp) {

    environment.init(envp, ".env");
    fileName = filesystem::path(__FILE__).filename();

    // TODO: implement print with less dataloss
    //       - method one: print everything 'manually'
    //       - method two: update only atom lines in the document
    printDocument(std::cout, "1a0s");
    return 0; // ignore remaining tests
    //
    //


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

std::string getPrefix(const std::string& key) {
    auto pos = key.find(".");
    if (pos == key.npos) {
        throw runtime_error(std::format("no prefix in '{}'", key));
    }
    return key.substr(0, pos);
}

std::string getSuffix(const std::string& key) {
    auto pos = key.find(".");
    if (pos == key.npos) {
        throw runtime_error(std::format("no prefix in '{}'", key));
    }
    return key.substr(pos + 1);
}

void printDocument(std::ostream& outputStream, const std::string& pdbCode) {

    auto inputPath = environment.get("PDB_CIF_DIR");
    inputPath += (string("/") + pdbCode[1] + pdbCode[2]) + "/" + pdbCode + "_updated.cif.gz";
    std::cerr << "reading: " << inputPath << " ..." << std::endl;

    auto document = gemmi::cif::read(gemmi::MaybeGzipped(inputPath));

    //
    // Update atom lines
    //

    auto oldBlock = document.blocks[0];
    // oldBlock.find_loop(); // col
    if (!oldBlock.has_mmcif_category("_atom_site")) {
        throw runtime_error("_atom_site not found");
    }
    auto atomLoopItem = *oldBlock.find_loop_item("_atom_site.id");
    auto atomTable = oldBlock.item_as_table(atomLoopItem);
    std::vector<std::string> columns;
    for (auto& tag : atomTable.tags()) {
        columns.emplace_back(getSuffix(tag.data()));
    }
    // auto& newLoop = document.blocks[0].init_mmcif_loop("_atom_site.", { "id", "name" });
    auto& newLoop = document.blocks[0].init_mmcif_loop("_atom_site.", columns);

    for (auto atom : atomTable) {
        // outputStream << atom.row_index << std::endl;
        std::vector<std::string> atomLine;
        for (auto& value : atom) {
            atomLine.emplace_back(value);
        }
        // TODO: replace this part by coords update
        // _atom_site.Cartn_x
        // _atom_site.Cartn_y
        // _atom_site.Cartn_z

        newLoop.add_row(atom);
    }


    // Printing...

    for (auto& block : document.blocks) {
        outputStream << std::format("data_{}", block.name) << std::endl;
        std::string lastPrefix{"data_"};
        for (auto& item : block.items) {
            // outputStream << item.line_number << std::endl;
            if (item.type == gemmi::cif::ItemType::Pair) {
                auto pair = item.pair;
                const std::string currentPrefix = getPrefix(pair[0]);
                if (currentPrefix != lastPrefix) {
                    lastPrefix = currentPrefix;
                    outputStream << "#" << std::endl;
                }
                const char separator = pair[1][0] == ';' ? '\n' : ' ';
                outputStream << std::format("{}{}{}", pair[0], separator, pair[1]) << std::endl;
            } else if (item.type == gemmi::cif::ItemType::Loop) {
                outputStream << "#" << std::endl;
                auto loop = item.loop;
                auto table = block.item_as_table(item);
                outputStream << "loop_" << std::endl;
                for (const auto& tag : table.tags()) {
                    outputStream << tag.data() << std::endl;
                }
                int colNum = table.tags().size();
                int column = 0;
                for (auto value : loop.values) {
                    column++;
                    outputStream << value;
                    if (column % colNum == 0) {
                        outputStream << std::endl;
                    } else {
                        outputStream << " ";
                    }
                }
                lastPrefix = "loop_";
            } else {
                outputStream << "unexpected type" << std::endl;
            }
        }
        outputStream << "#" << std::endl;
    }

}
