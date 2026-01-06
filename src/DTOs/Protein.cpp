// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

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


        void printDocument(std::ostream& outputStream, Tmdet::VOs::Protein& protein);

        std::stringstream sstream;
        printDocument(sstream, protein);

        if (path.ends_with(".gz")) {
            Tmdet::Helpers::Gzip::writeFile(path, sstream.str());
        } else {
            std::ofstream outCif(path);
            outCif << sstream.str();
        }

    }

    void printDocument(std::ostream& outputStream, Tmdet::VOs::Protein& protein) {

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
                    const std::string currentPrefix = Tmdet::Utils::CifUtil::getPrefix(pair[0]);
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
                    auto loop = item.loop;
                    auto table = block.item_as_table(item);
                    if ((int)table.length() == 0) {
                        // skip empty loops, chimera parser fails on them
                        continue;
                    }
                    outputStream << "#" << std::endl;
                    outputStream << "loop_" << std::endl;
                    // process the loop as table
                    for (const auto& tag : table.tags()) {
                        outputStream << tag.data() << std::endl;
                    }
                    // number of columns in the table
                    int colNum = table.tags().size();
                    int column = 0;
                    // previousWasLongText flag to avoid double std::endl after long value
                    bool previousWasLongText = false;
                    // iterate on values in the loop object;
                    // 'column' helps to identify end of rows
                    for (auto value : loop.values) {
                        column++;
                        // print value of a column
                        if (value[0] == ';' && !previousWasLongText) {
                            // insert new line before long string separator,
                            // except the previous long value already printed it
                            outputStream << std::endl;
                            previousWasLongText = true;
                        } else {
                            previousWasLongText = false;
                        }
                        outputStream << value;
                        // if new row begins or last char is ; (so it is a long value)
                        if (column % colNum == 0 || *(value.end() - 1) == ';') {
                            outputStream << std::endl;
                        } else {
                            // separator between column values
                            outputStream << " ";
                        }
                    }
                    lastPrefix = "loop_";
                } else if (item.type == gemmi::cif::ItemType::Erased) {
                    // NOTE: only pair expected as erased item (see addMembraneEntity function)
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

    Tmdet::VOs::Protein Protein::get(const std::string& inputPath, const int modelIndex) {
        Tmdet::VOs::Protein protein;
        protein.getStructure(inputPath);
        protein.code = protein.gemmi.name;
        protein.inputFile = inputPath;
        Tmdet::Helpers::String::toLower(protein.code);
        protein.modelIndex = (modelIndex>=(int)protein.gemmi.models.size()?0:modelIndex);
        remove_hydrogens(protein.gemmi.models[protein.modelIndex]);
        remove_ligands_and_waters(protein.gemmi.models[protein.modelIndex]);
        remove_alternative_conformations(protein.gemmi.models[protein.modelIndex]);
        //protein.gemmi.models.resize(1);

        int chainIdx = 0;
        for(auto& chain: protein.gemmi.models[protein.modelIndex].chains) {
            protein.chains.emplace_back(Tmdet::DTOs::Chain::get(protein.gemmi,chain,chainIdx));
            chainIdx++;
        }
        if (protein.gemmi.cell.a == 0.0 && protein.gemmi.cell.b == 0.0 && protein.gemmi.cell.c == 0.0 ) {
            protein.gemmi.cell.a = 1.0;
            protein.gemmi.cell.b = 1.0;
            protein.gemmi.cell.c = 1.0;
        }
        protein.neighbors = gemmi::NeighborSearch(protein.gemmi.models[protein.modelIndex], protein.gemmi.cell, 9);
        protein.neighbors.populate();
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
                protein.bioMatrix.deletedChainIds.push_back(chainId);
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
            ret += Tmdet::DTOs::SecStrVec::toString(protein,secStrVec);
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
                addCurvedMembraneAtoms(protein, membrane, ret);
            }
        }
        return ret;
    }

    void Protein::addPlaneMembraneAtoms(Tmdet::VOs::Protein& protein, const Tmdet::VOs::Membrane& membrane, std::vector<gemmi::Vec3>& ret) {
        double R = membrane.membraneRadius + 10;
        int sign = 1;
        for (double x=-R; x<=R; x+=6) {
            for (double y=-R + sign * 1.5; y<=R + sign * 1.5; y+=6) {
                if (sqrt(x*x+y*y) < R) {
                    ret.push_back(gemmi::Vec3(x,y,membrane.origo+membrane.halfThickness));
                    ret.push_back(gemmi::Vec3(x,y,membrane.origo-membrane.halfThickness));
                }
                sign *= -1;
            }
        }
    }

    void Protein::addCurvedMembraneAtoms(Tmdet::VOs::Protein& protein, const Tmdet::VOs::Membrane& membrane, std::vector<gemmi::Vec3>& ret) {
        double R = membrane.membraneRadius;
        int sign = 1;
        for (double x=-R; x<=R; x+=6) {
            for (double y=-R + sign * 1.5; y<=R + sign * 1.5; y+=6) {
                if (sqrt(x*x+y*y) < R) {
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
                sign *= -1;
            }
        }
    }
}
