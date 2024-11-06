#include <sstream>
#include <string>
#include <gemmi/modify.hpp>
#include <gemmi/polyheur.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/to_mmcif.hpp>
#include <Config.hpp>
#include <DTOs/Chain.hpp>
#include <DTOs/SecStrVec.hpp>
#include <DTOs/Protein.hpp>
#include <Helpers/Gzip.hpp>
#include <Helpers/String.hpp>
#include <System/FilePaths.hpp>
#include <System/Logger.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Chain.hpp>
#include <VOs/Residue.hpp>
#include <VOs/Atom.hpp>

namespace Tmdet::DTOs {

    void Protein::writeCif(Tmdet::VOs::Protein& protein, const std::string& path) {

        logger.debug("Processing Protein::writeCif()");

        void printDocument(std::ostream& outputStream, Tmdet::VOs::Protein& protein);

        std::stringstream sstream;
        printDocument(sstream, protein);

        if (path.ends_with(".gz")) {
            Tmdet::Helpers::Gzip::writeFile(path, sstream.str());
        } else {
            std::ofstream outCif(path);
            outCif << sstream.str();
        }

        logger.debug(" Document is written into {}", path);
        logger.debug(" Processed Protein::writeCif()");
    }

    void printDocument(std::ostream& outputStream, Tmdet::VOs::Protein& protein) {

        //
        // Helper functions
        //

        auto getPrefix = [](const std::string& key) -> std::string {
            auto pos = key.find(".");
            if (pos == key.npos) {
                throw std::runtime_error(std::format("no prefix in '{}'", key));
            }
            return key.substr(0, pos);
        };

        auto getSuffix = [](const std::string& key) -> std::string {
            auto pos = key.find(".");
            if (pos == key.npos) {
                throw std::runtime_error(std::format("no prefix in '{}'", key));
            }
            return key.substr(pos + 1);
        };

        /////////////////////////////////////////////////////////////
        //
        // Main logic
        //
        /////////////////////////////////////////////////////////////


        //
        // Update atom lines in the document
        //
        auto& document = protein.document;
        auto oldBlock = document.blocks[0];

        if (!oldBlock.has_mmcif_category("_atom_site")) {
            throw std::runtime_error("_atom_site not found");
        }
        auto atomLoopItem = *oldBlock.find_loop_item("_atom_site.id");
        auto atomTable = oldBlock.item_as_table(atomLoopItem);
        // collect column names and store indecies of x,y,z coord columns
        std::vector<std::string> columns;
        int xIndex = -1;
        int yIndex = -1;
        int zIndex = -1;

        // set column indecies
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
            }
            // } else if (colName == "id") {
            //     serialIndex = colIndex;
            // }
            colIndex++;
        }
        auto& newLoop = document.blocks[0].init_mmcif_loop("_atom_site.", columns);

        // update atom coords
        auto updateCoords = [protein](std::vector<std::string>& values, int xColumn, int yColumn, int zColumn) {

            double x = std::stod(values[xColumn]);
            double y = std::stod(values[yColumn]);
            double z = std::stod(values[zColumn]);

            gemmi::Vec3 pos{x, y, z};
            pos = protein.tmatrix.rot.multiply(pos);

            x = protein.tmatrix.trans.x + pos.x;
            y = protein.tmatrix.trans.y + pos.y;
            z = protein.tmatrix.trans.z + pos.z;

            values[xColumn] = std::format("{:.3f}", x);
            values[yColumn] = std::format("{:.3f}", y);
            values[zColumn] = std::format("{:.3f}", z);
        };
        for (auto atom : atomTable) {
            // outputStream << atom.row_index << std::endl;
            std::vector<std::string> atomLine;
            for (auto& value : atom) {
                atomLine.emplace_back(value);
            }
            // transform atom position
            updateCoords(atomLine, xIndex, yIndex, zIndex);
            // overwrite atom position id in the Table
            atom[xIndex] = atomLine[xIndex];
            atom[yIndex] = atomLine[yIndex];
            atom[zIndex] = atomLine[zIndex];
            // add atom row to the loop
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

    Tmdet::VOs::Protein Protein::get(const std::string& inputPath) {
        DEBUG_LOG("Processing Protein::get()");
        Tmdet::VOs::Protein protein;
        protein.getStructure(inputPath);
        protein.code = protein.gemmi.name;
        Tmdet::Helpers::String::toLower(protein.code);
        remove_hydrogens(protein.gemmi.models[0]);
        remove_ligands_and_waters(protein.gemmi.models[0]);
        remove_alternative_conformations(protein.gemmi.models[0]);
        protein.gemmi.models.resize(1);

        int chainIdx = 0;
        for(auto& chain: protein.gemmi.models[0].chains) {
            protein.chains.emplace_back(Tmdet::DTOs::Chain::get(protein.gemmi,chain,chainIdx));
            chainIdx++;
        }
        protein.neighbors = gemmi::NeighborSearch(protein.gemmi.models[0], protein.gemmi.cell, 9);
        protein.neighbors.populate();
        DEBUG_LOG(" Processed Protein::get()");
        return protein;
    }

    void Protein::unselectAntiBodyChains(Tmdet::VOs::Protein& protein) {
        for (auto& chain : protein.chains) {
            if (!protein.polymerNames.contains(chain.entityId)) {
                continue;
            }
            auto name = protein.polymerNames[chain.entityId];
            for (auto& filter : Tmdet::ANTIBODY_NAMES) {
                if (Tmdet::Helpers::String::toUpper(name).find(filter) != name.npos) {
                    INFO_LOG("Unselecting Ab chain: {}",chain.id);
                    chain.selected = false;
                    break;
                }
            }
        }
    }

    void Protein::unselectChains(const std::string& chainIds, Tmdet::VOs::Protein& protein) {
        for(auto chainId: Tmdet::Helpers::String::explode(",",chainIds)) {
            if (int chainIdx = protein.searchChainById(chainId); chainIdx != -1) {
                protein.chains[chainIdx].selected = false;
                INFO_LOG("Unselecting chain: {}",chainId);
            }
            else {
                WARN_LOG("Could not find chain: {}",chainId);
            }
        }
    }

    std::string Protein::toString(const Tmdet::VOs::Protein& protein) {
        std::string ret = "";
        for(const auto& chain: protein.chains) {
            ret += Tmdet::DTOs::Chain::toString(chain);
        }
        for (const auto& secStrVec: protein.secStrVecs) {
            ret += Tmdet::DTOs::SecStrVec::toString(secStrVec);
        }
        return ret;
    }

    std::vector<gemmi::Vec3> Protein::addMembraneAtoms(Tmdet::VOs::Protein& protein) {
        std::vector<gemmi::Vec3> ret;
        for (const auto& membrane: protein.membranes) {
            if (membrane.type.isPlane()) {
                addPlaneMembraneAtoms(protein, membrane, ret);
            }
            else {
                addBlendedMembraneAtoms(protein, membrane, ret);
            }
        }
        return ret;
    }

    void Protein::addPlaneMembraneAtoms(Tmdet::VOs::Protein& protein, const Tmdet::VOs::Membrane& membrane, std::vector<gemmi::Vec3>& ret) {
        for (double x=-membrane.membraneRadius; x<=membrane.membraneRadius; x+=3) {
            for (double y=-membrane.membraneRadius; y<=membrane.membraneRadius; y+=3) {
                if (sqrt(x*x+y*y) < membrane.membraneRadius) {
                    ret.push_back(gemmi::Vec3(x,y,membrane.halfThickness));
                    ret.push_back(gemmi::Vec3(x,y,-membrane.halfThickness));
                }
            }
        }
    }

    void Protein::addBlendedMembraneAtoms(Tmdet::VOs::Protein& protein, const Tmdet::VOs::Membrane& membrane, std::vector<gemmi::Vec3>& ret) {
        for (double x=-membrane.membraneRadius; x<=membrane.membraneRadius; x+=3) {
            for (double y=-membrane.membraneRadius; y<=membrane.membraneRadius; y+=3) {
                if (sqrt(x*x+y*y) < membrane.membraneRadius) {
                    double r = membrane.sphereRadius + membrane.halfThickness;
                    double z = sqrt(r*r - x*x - y*y);
                    ret.push_back(gemmi::Vec3(x,y,z));
                    r = membrane.sphereRadius - membrane.halfThickness;
                    z = sqrt(r*r - x*x - y*y);
                    ret.push_back(gemmi::Vec3(x,y,z));
                }
            }
        }
    }
}
