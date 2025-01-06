// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <format>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <DTOs/Protein.hpp>
#include <Utils/CifUtil.hpp>

using namespace std;

namespace Tmdet::Utils {

    std::string CifUtil::getSuffix(const std::string& tag) {
        auto pos = tag.find(".");
        if (pos == tag.npos) {
            throw std::runtime_error(std::format("no prefix in '{}'", tag));
        }
        return tag.substr(pos + 1);
    }

    std::string addMembraneEntity(gemmi::cif::Block& block) {
        if (!block.has_mmcif_category("_entity")) {
            throw std::runtime_error("_entity not found");
        }

        if (block.has_mmcif_category("_pdbx_struct_assembly_gen")) {
            // TODO: what if this category is a loop?
            auto asymIdListPair = block.find_pair("_pdbx_struct_assembly_gen.asym_id_list");
            if (asymIdListPair != nullptr) {
                std::string newList = std::format("{},{}", (*asymIdListPair)[1], CifUtil::TMDET_MEMBRANE_ASYM_ID);
                block.set_pair("_pdbx_struct_assembly_gen.asym_id_list", newList);
            }
        }

        // TODO: find_mmcif_category ???
        auto entityItemPtr = block.find_loop_item("_entity.id");
        if (entityItemPtr == nullptr) {
            entityItemPtr = block.find_pair_item("_entity.id");
        }
        auto entityItem = *entityItemPtr;

        auto entityTable = entityItem.type == gemmi::cif::ItemType::Loop
            ? block.item_as_table(entityItem)
            // TODO: use each columns
            : block.find("_entity.", { "id", "type", "pdbx_description" });


        std::vector<std::string> columns{
            "id",
            "type",
            "pdbx_description",
        };

        // Collect row values
        std::vector<std::vector<std::string>> newRows;
        for (const auto& row : entityTable) {
            newRows.push_back({ row.at(0), row.at(1), row.at(2) });
        }

        // Overwrite _entity category by a new loop
        auto& newLoop = block.init_mmcif_loop("_entity.", columns);
        for (const auto& row : newRows) {
            newLoop.add_row(row);
        }

        std::string newEntityId = std::format("{:d}", newRows.size() + 1);
        std::vector<std::string> values{
            newEntityId,
            ".",
            gemmi::cif::quote("TMDET MEMBRANE REPRESENTATION"),
        };
        newLoop.add_row(values);

        return newEntityId;
    }

    void CifUtil::prepareDocumentBlock(Tmdet::VOs::Protein& protein) {

        auto& document = protein.document;
        auto& oldBlock = document.blocks[0];

        auto entityId = addMembraneEntity(oldBlock);

        if (!oldBlock.has_mmcif_category("_atom_site")) {
            throw std::runtime_error("_atom_site not found");
        }
        auto atomLoopItem = *oldBlock.find_loop_item("_atom_site.id");
        auto atomTable = oldBlock.item_as_table(atomLoopItem);
        // collect column names and store indecies of x,y,z coord columns
        std::vector<std::string> columns;
        columns.clear();
        int xIndex = -1;
        int yIndex = -1;
        int zIndex = -1;
        int serialIndex = -1;
        int typeSymbolIndex = -1;
        int labelCompIdIndex = -1;
        int labelAsymIdIndex = -1;
        int labelEntityIdIndex = -1;
        int labelSeqIdIndex = -1;
        int occupancyIndex = -1;
        int bIsoIndex = -1;
        int pdbxPDBModelNumIndex = -1;

        // set column indecies
        int colIndex = 0;
        for (auto& tag : atomTable.tags()) {
            auto colName = Tmdet::Utils::CifUtil::getSuffix(tag.data());
            columns.emplace_back(colName);
            if (colName == "Cartn_x") {
                xIndex = colIndex;
            } else if (colName == "Cartn_y") {
                yIndex = colIndex;
            } else if (colName == "Cartn_z") {
                zIndex = colIndex;
            } else if (colName == "id") {
                serialIndex = colIndex;
            } else if (colName == "type_symbol") {
                typeSymbolIndex = colIndex;
            } else if (colName == "label_comp_id") {
                labelCompIdIndex = colIndex;
            } else if (colName == "label_asym_id") {
                labelAsymIdIndex = colIndex;
            } else if (colName == "label_entity_id") {
                labelEntityIdIndex = colIndex;
            } else if (colName == "label_seq_id") {
                labelSeqIdIndex = colIndex;
            } else if (colName == "occupancy") {
                occupancyIndex = colIndex;
            } else if (colName == "B_iso_or_equiv") {
                bIsoIndex = colIndex;
            } else if (colName == "pdbx_PDB_model_num") {
                pdbxPDBModelNumIndex = colIndex;
            }

            colIndex++;
        }
        // Init new loop for ATOM/HETATM lines
        auto& newLoop = document.blocks[0].init_mmcif_loop("_atom_site.", columns);

        // lambda function to update atom coords - including transformation
        auto updateCoords = [&](std::vector<std::string>& values, int xColumn, int yColumn, int zColumn) {

            double x = std::stod(values[xColumn]);
            double y = std::stod(values[yColumn]);
            double z = std::stod(values[zColumn]);

            gemmi::Vec3 pos{x, y, z};
            protein.tmatrix.transform(pos);
            values[xColumn] = std::format("{:.3f}", pos.x);
            values[yColumn] = std::format("{:.3f}", pos.y);
            values[zColumn] = std::format("{:.3f}", pos.z);
        };

        // Iterate on existing atoms
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

        // Add atoms for membrane representation
        // TODO:

        // colIndex == number of columns
        const auto& membraneAtoms = Tmdet::DTOs::Protein::addMembraneAtoms(protein);

        int nextSerialId = atomTable.length() + 1;
        for (const auto& atomPosition : membraneAtoms) {
            std::vector<std::string> atomLine(colIndex, ".");
            atomLine[0] = "HETATM";
            atomLine[serialIndex] = std::format("{:d}", nextSerialId);
            atomLine[typeSymbolIndex] = "AG";
            atomLine[labelCompIdIndex] = "AG";
            atomLine[labelAsymIdIndex] = TMDET_MEMBRANE_ASYM_ID;
            atomLine[labelEntityIdIndex] = entityId;
            atomLine[labelSeqIdIndex] = "1"; // This is a single giant residue
            atomLine[occupancyIndex] = "1.00";
            atomLine[bIsoIndex] = "0.00";
            atomLine[pdbxPDBModelNumIndex] = "1";
            atomLine[xIndex] = std::format("{:.3f}", atomPosition.x);
            atomLine[yIndex] = std::format("{:.3f}", atomPosition.y);
            atomLine[zIndex] = std::format("{:.3f}", atomPosition.z);

            newLoop.add_row(atomLine);
            nextSerialId++;
        }

    }
}
