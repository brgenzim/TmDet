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
#include <Utils/CifUtil.hpp>
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

        /////////////////////////////////////////////////////////////
        //
        // Main logic
        //
        /////////////////////////////////////////////////////////////

        //
        // Update _atom_site loop, _entity and _pdbx_struct_assembly_gen categories
        //
        Tmdet::Utils::CifUtil::prepareDocumentBlock(protein);

        // Printing...
        auto& document = protein.document;

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
                        if (value[0] == ';') {
                            // insert new line before long string separator
                            outputStream << std::endl;
                        }
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
                } else if (item.type == gemmi::cif::ItemType::Erased) {
                    // NOTE: only pair expected as erased item (see addMembraneEntity function)
                    // outputStream << std::format("# TMDET warning: erased pair: {} {}", item.pair[0], item.pair[1]) << std::endl;
                    // DEBUG_LOG("# TMDET warning: erased pair: {} {}", item.pair[0], item.pair[1]);
                    continue;
                } else {
                    std::unordered_map<gemmi::cif::ItemType, std::string> types{
                        { gemmi::cif::ItemType::Frame, "Frame" },
                        { gemmi::cif::ItemType::Comment, "Comment" },
                    };

                    outputStream << "# TMDET warning: unexpected type: " << types[item.type] << std::endl;
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
        for (double x=-membrane.membraneRadius; x<=membrane.membraneRadius; x+=5) {
            for (double y=-membrane.membraneRadius; y<=membrane.membraneRadius; y+=5) {
                if (sqrt(x*x+y*y) < membrane.membraneRadius) {
                    ret.push_back(gemmi::Vec3(x,y,membrane.halfThickness));
                    ret.push_back(gemmi::Vec3(x,y,-membrane.halfThickness));
                }
            }
        }
    }

    void Protein::addBlendedMembraneAtoms(Tmdet::VOs::Protein& protein, const Tmdet::VOs::Membrane& membrane, std::vector<gemmi::Vec3>& ret) {
        for (double x=-membrane.membraneRadius; x<=membrane.membraneRadius; x+=5) {
            for (double y=-membrane.membraneRadius; y<=membrane.membraneRadius; y+=5) {
                if (sqrt(x*x+y*y) < membrane.membraneRadius) {
                    double r = membrane.sphereRadius + membrane.halfThickness;
                    if (r*r > x*x + y*y) {
                        double z = sqrt(r*r - x*x - y*y) + membrane.origo;
                        ret.push_back(gemmi::Vec3(x,y,z));
                    }
                    r = membrane.sphereRadius - membrane.halfThickness;
                    if (r*r > x*x + y*y) {
                        double z = sqrt(r*r - x*x - y*y) + membrane.origo;
                        ret.push_back(gemmi::Vec3(x,y,z));
                    }
                }
            }
        }
    }
}
