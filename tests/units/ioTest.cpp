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
void printDocument(std::ostream& outputStream, Tmdet::ValueObjects::Protein& protein);

string fileName;
Tmdet::System::Environment environment;
Tmdet::System::Logger logger;

int main(int argc, char *argv[], char **envp) {

    environment.init(envp, ".env");
    fileName = filesystem::path(__FILE__).filename();

    std::ofstream logFile("build/ioTest.log");
    logger.addStream(logFile);
    logger.setLevel(Tmdet::System::level::debug);

    // TODO: implement print with less dataloss
    //       [x] method one: print everything 'manually'
    //       - method two: update only atom lines in the document





    // Init
    auto protein = createTmdetStruct("1a0s");

    // Test case 1 - CIF write test
    {
        protein.tmatrix.rot = {
            1, 0, 0,
            0, 0, -1,
            0, 1, 0
        };
        Tmdet::DTOs::Protein::transform(protein);
        printDocument(logFile, protein);
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

void printDocument(std::ostream& outputStream, Tmdet::ValueObjects::Protein& protein) {

    //
    // Helper functions
    //

    auto getPrefix = [](const std::string& key) -> std::string {
        auto pos = key.find(".");
        if (pos == key.npos) {
            throw runtime_error(std::format("no prefix in '{}'", key));
        }
        return key.substr(0, pos);
    };

    auto getSuffix = [](const std::string& key) -> std::string {
        auto pos = key.find(".");
        if (pos == key.npos) {
            throw runtime_error(std::format("no prefix in '{}'", key));
        }
        return key.substr(pos + 1);
    };

    /////////////////////////////////////////////////////////////
    //
    // Main logic
    //
    /////////////////////////////////////////////////////////////

    auto& document = protein.document;

    std::vector<gemmi::Position> atomPositions;
    int nextSerial = 1;
    for (auto& chain: protein.chains) {
        for (auto& residue: chain.residues) {
            for (auto& atom: residue.atoms) {
                if (atom.gemmi.serial != nextSerial) {
                    throw runtime_error(std::format(
                        "Unexpected atom serial number: {} / expected number: {}", atom.gemmi.serial, nextSerial
                    ));
                }
                nextSerial++;
            }
        }
    }


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
    // collect column names and store indecies of x,y,z coord columns
    std::vector<std::string> columns;
    int xIndex = -1;
    int yIndex = -1;
    int zIndex = -1;
    int serialIndex = -1;

    int colIndex = 0;
    for (auto& tag : atomTable.tags()) {
        auto colName = getSuffix(tag.data());
        columns.emplace_back(colName);
        if (colName == "Cartn_x") {
            xIndex = colIndex;
        } else if (colName == "Cartn_y") {
            yIndex = colIndex;
        } else if (colName == "Cartn_z") {
            zIndex = colIndex;
        } else if (colName == "id") {
            serialIndex = colIndex;
        }
        colIndex++;
    }
    auto& newLoop = document.blocks[0].init_mmcif_loop("_atom_site.", columns);

    // set column indecies

    for (auto atom : atomTable) {
        // outputStream << atom.row_index << std::endl;
        std::vector<std::string> atomLine;
        for (auto& value : atom) {
            atomLine.emplace_back(value);
        }
        // TODO: replace this part by real transformation
        auto serialId = std::stoi(atom[serialIndex]);
        const auto& position = atomPositions[serialId - 1];
        atom[xIndex] = std::format("{:.3f}", position.x);
        atom[yIndex] = std::format("{:.3f}", position.y);
        atom[zIndex] = std::format("{:.3f}", position.z);

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
                // get category prefix
                const std::string currentPrefix = getPrefix(pair[0]);
                if (currentPrefix != lastPrefix) {
                    // if prefix changes print category delimiter
                    lastPrefix = currentPrefix;
                    outputStream << "#" << std::endl;
                }
                // if long string value (with ';' - boundaries)
                // then print '\n' between tag and its value
                // else just use a space
                const char separator = pair[1][0] == ';' ? '\n' : ' ';
                outputStream << std::format("{}{}{}", pair[0], separator, pair[1]) << std::endl;
            } else if (item.type == gemmi::cif::ItemType::Loop) {
                outputStream << "#" << std::endl;
                auto loop = item.loop;
                auto table = block.item_as_table(item);
                outputStream << "loop_" << std::endl;
                // process the loop as table
                for (const auto& tag : table.tags()) {
                    outputStream << tag.data() << std::endl;
                }
                // number of columns in the table
                int colNum = table.tags().size();
                int column = 0;
                // iterate on values in the loop object;
                // 'column' helps to identify end of rows
                for (auto value : loop.values) {
                    column++;
                    // print value of a column
                    outputStream << value;
                    // if new row begins
                    if (column % colNum == 0) {
                        outputStream << std::endl;
                    } else {
                        // separator between column values
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
