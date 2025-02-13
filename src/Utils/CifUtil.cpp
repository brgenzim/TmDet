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
#include <Services/ChemicalComponentDirectoryService.hpp>

using namespace std;

namespace Tmdet::Utils {

    std::string CifUtil::getSuffix(const std::string& tag) {
        auto pos = tag.find(".");
        if (pos == tag.npos) {
            throw std::runtime_error(std::format("no prefix in '{}'", tag));
        }
        return tag.substr(pos + 1);
    }

    void updateStructAssemblyGen(gemmi::cif::Block& block) {
        // local function
        auto appendMembraneId = [](const std::string& idList) -> std::string {
            std::string result{idList};
            if (idList[0] == ';') {
                auto insertionPosition = (idList.ends_with("\r\n;"))
                    ? idList.size() - 3
                    // otherwise it ends with "\n;"
                    : idList.size() - 2;
                result.insert(insertionPosition, "," + CifUtil::TMDET_MEMBRANE_ASYM_ID);
            } else {
                result += ("," + CifUtil::TMDET_MEMBRANE_ASYM_ID);
            }
            return result;
        };

        // relevant columns (others will be excluded from the new loop):
        std::vector<std::string> columns{ "assembly_id", "oper_expression", "asym_id_list" };
        auto asymGenTable = block.find("_pdbx_struct_assembly_gen.", columns);

        if (asymGenTable.size() == 0) {
            return;
        }

        // Collect row values
        std::vector<std::vector<std::string>> newRows;
        for (const auto& row : asymGenTable) {
            newRows.push_back({ row.at(0), row.at(1), row.at(2) });
        }
        // _pdbx_struct_assembly_gen.assembly_id
        // _pdbx_struct_assembly_gen.oper_expression
        // _pdbx_struct_assembly_gen.asym_id_list
        auto& newGenLoop = block.init_mmcif_loop("_pdbx_struct_assembly_gen.", columns);
        for (auto& row : newRows) {
            // add TM_ entity to each assembly
            row[2] = appendMembraneId(row[2]);
            newGenLoop.add_row(row);
        }
    }

    std::string addMembraneEntity(gemmi::cif::Block& block) {
        if (!block.has_mmcif_category("_entity")) {
            throw std::runtime_error("_entity not found");
        }

        if (block.has_mmcif_category("_pdbx_struct_assembly_gen")) {
            updateStructAssemblyGen(block);
        }

        // get entity list as table

        // relevant columns (other columns will be omitted in the output)
        std::vector<std::string> columns{
            "id",
            "type",
            "pdbx_description",
        };
        std::vector<std::vector<std::string>> newRows;

        {
            auto entityTable = block.find("_entity.", columns);

            // Collect row values
            for (const auto& row : entityTable) {
                newRows.push_back({ row.at(0), row.at(1), row.at(2) });
            }

            // loop with these columns not found
            if (newRows.size() == 0) {
                // drop "pdbx_description";
                // Probably tmdet input is a CIF file
                // converted by 'gemmi convert' from an ENT file.
                columns.erase(columns.end());

                auto entityTable = block.find("_entity.", columns);
                // Collect row values and append unknown pdbx_description
                for (const auto& row : entityTable) {
                    newRows.push_back({ row.at(0), row.at(1), "." });
                }
                // We will use these unknown descriptions;
                // so we can reconstruct original entities
                // from the input CIF in the overwrite block.
                columns.push_back("pdbx_description");
            }
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

    /**
     * @brief gemmi convert from pdb to cif does not add type details to
     *        _chem_comp category. Therefore Mol* cannot apply cartoon
     *        representation on polymer chains and displays everything with
     *        space-ball representation.
     *        This function update the given document block by the missing
     *        details if the first type value is just a single "." (dot).
     */
    void updateChemicalComponentLoopIfNeeded(gemmi::cif::Block& block) {

        // check whether this document resulted by pdb-cif conversion
        auto chemCompTable = block.find("_chem_comp.", { "id", "type" });
        std::vector<std::string> componentIds;
        {
            bool isTypeMissing = false;
            for (const auto& row : chemCompTable) {
                componentIds.push_back(row.at(0));
                if (!isTypeMissing && row.at(1) == ".") {
                    isTypeMissing = true;
                }
            }
            if (!isTypeMissing) {
                // _chemp_comp category has type information, update isn't needed
                return;
            }
        }

        //
        // correction of _chem_comp
        //

        std::vector<std::string> columns{ "id", "type", "name", "one_letter_code",
            "three_letter_code", "formula", "formula_weight" };
        auto& newChemLoop = block.init_mmcif_loop("_chem_comp.", columns);

        for (const auto& threeLetterCode : componentIds) {
            const auto& row = Tmdet::Services::ChemicalComponentDirectoryService::getChemicalComponentInfo(threeLetterCode, columns);
            newChemLoop.add_row(row);
        }
    }

    void CifUtil::prepareDocumentBlock(Tmdet::VOs::Protein& protein) {

        int modelIndex = protein.modelIndex;

        auto& document = protein.document;
        auto& oldBlock = document.blocks[0];

        auto entityId = addMembraneEntity(oldBlock);
        updateChemicalComponentLoopIfNeeded(oldBlock);

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
        const std::string modelNumber = std::to_string(protein.gemmi.models[modelIndex].num);
        int nextSerialId = 1;
        for (auto atom : atomTable) {

            // skip unneeded models
            if (atom[pdbxPDBModelNumIndex] != modelNumber) {
                continue;
            }

            atom[serialIndex] = std::format("{:d}", nextSerialId);

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
            nextSerialId++;
        }

        // Add atoms for membrane representation

        // colIndex == number of columns
        const auto& membraneAtoms = Tmdet::DTOs::Protein::addMembraneAtoms(protein);

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
            atomLine[pdbxPDBModelNumIndex] = modelNumber;
            atomLine[xIndex] = std::format("{:.3f}", atomPosition.x);
            atomLine[yIndex] = std::format("{:.3f}", atomPosition.y);
            atomLine[zIndex] = std::format("{:.3f}", atomPosition.z);

            newLoop.add_row(atomLine);
            nextSerialId++;
        }

    }

}
