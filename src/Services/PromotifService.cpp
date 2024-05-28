#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cstdio>
#include <memory>

#include <gemmi/cif.hpp>
#include <gemmi/cifdoc.hpp>
#include <gemmi/mmcif.hpp>

#include <Services/PromotifService.hpp>

#define READ_BUFFER_SIZE 200
#define CMD_TEMPLATE "cat {} | process_entry -input /dev/stdin -input_format cif -output /dev/stdout -output_format cif -keep_original_numbering"
#define SECONDARY_STRUCTURE_COLUMN 17

using namespace std;

namespace Tmdet::Services::PromotifService {

    struct SecondaryStruct {
        char code;
        string labelAsymId;
        int sequenceBegin;
        int sequenceEnd;
    };

    static std::map<std::string, std::string> parseOutput(const std::string& promotifOutput);
    static std::string exec(const std::string& cmd);
    static std::vector<string> getChainSequence(const gemmi::Structure& structure, const string& chainName);
    static std::map<string, std::vector<SecondaryStruct>> parseStructConf(gemmi::cif::Block& block);
    static std::map<string, std::vector<SecondaryStruct>> parseSheetRange(gemmi::cif::Block& block);
    static void updateDsspStrings(std::map<string, string>& dsspStrings, std::map<string, std::vector<SecondaryStruct>>& secondaryStructs);
    static char getSecondaryStructureType(const string& mmCifType);


    std::map<std::string, std::string> process(const std::string& cifPath) {

        // TODO: get RCSBROOT env var
        // TODO: check process_entry is in PATH env var

        std::string cmd = std::string(CMD_TEMPLATE);
        // NOTE std::format is not supported in gcc 11 (Ubuntu 22.04)
        if (cifPath.ends_with("gz")) {
            // use zcat instead of 'cat' defaulted in the command template
            cmd = "z" + cmd;
        }
        cmd = cmd.replace(cmd.find("{}"), 2, cifPath);
        std::string promotifOutput = exec(cmd);

        return parseOutput(promotifOutput);
    }

    static std::string exec(const std::string& cmd) {
        std::array<char, READ_BUFFER_SIZE> buffer;
        std::string result;
        std::shared_ptr<FILE> pipe(popen(cmd.c_str(), "r"), pclose);
        if (!pipe) {
            throw std::runtime_error("popen() failed! Command: " + cmd);
        }
        while (!feof(pipe.get())) {
            if (fgets(buffer.data(), READ_BUFFER_SIZE, pipe.get()) != nullptr)
                result += buffer.data();
        }
        return result;
    }

    static std::map<std::string, std::string> parseOutput(const std::string& promotifOutput) {

        std::map<string, string> dsspStrings;
        std::istringstream stream(promotifOutput);

        gemmi::cif::Document document = gemmi::cif::read_string(promotifOutput);
        gemmi::Structure structure = gemmi::make_structure(std::move(document));
        gemmi::cif::Block& block = document.blocks[0];

        auto structConfs = parseStructConf(block);
        auto sheetRanges = parseSheetRange(block);

        // init dssp strings
        for (const auto& chainId : structConfs) {
            const auto& sequence = getChainSequence(structure, chainId.first);
            string dssp(sequence.size(), '-');
            dsspStrings[chainId.first] =  dssp;
        }
        for (const auto& chainId : sheetRanges) {
            if (dsspStrings.contains(chainId.first)) {
                continue;
            }
            const auto& sequence = getChainSequence(structure, chainId.first);
            string dssp(sequence.size(), '-');
            dsspStrings[chainId.first] =  dssp;
        }

        updateDsspStrings(dsspStrings, structConfs);
        updateDsspStrings(dsspStrings, sheetRanges);

        for (auto& [key, value] : dsspStrings) {
            cout << key << ": " << value << endl;
        }
        return dsspStrings;
    }

    std::vector<string> getChainSequence(const gemmi::Structure& structure, const string& chainName) {

        std::vector<string> sequence;
        const gemmi::Chain* chain = structure.first_model().find_chain(chainName);
        if (chain == nullptr) {
            throw runtime_error("Chain '" + chainName + "' not found");
        }
        auto entityId = chain->residues[0].entity_id;
        for (const auto& entity : structure.entities) {
            if (entity.entity_type == gemmi::EntityType::Polymer
                && entity.name == entityId) {

                sequence = entity.full_sequence;
                break;
            }
        }
        return sequence;
    }

    char getSecondaryStructureType(const string& mmCifType) {
        char cType = '-';
        if (mmCifType == "BEND") {
            cType = 'S';
        } else if (mmCifType == "HELX_RH_3T_P") {
            cType = 'G';
        } else if (mmCifType == "HELX_RH_PI_P") {
            cType = 'I';
        } else if (mmCifType == "HELX_LH_PP_P") {
            cType = 'P';
        } else if (mmCifType == "HELX_RH_AL_P") {
            cType = 'H';
        // other helix types
        } else if (mmCifType.find("HELX") == 0) {
            cType = 'H';
        } else if (mmCifType == "OTHER") {
            cType = ' ';
        } else if (mmCifType == "STRN") {
            cType = 'B';
        } else if (mmCifType.find("TURN") == 0) {
            cType = 'T';
        } else {
            throw std::runtime_error("Unexpected mmCIF secondary strcuture type: '" + mmCifType + "'");
        }
        return cType;
    }

    std::map<string, std::vector<SecondaryStruct>> parseStructConf(gemmi::cif::Block& block) {
        std::map<string, std::vector<SecondaryStruct>> secondaryStructures;
        // TODO: check: has_mmcif_category
        auto structConfLoop = block.find_loop_item("_struct_conf.id")->loop;
        int loopLength = structConfLoop.length();
        // int idCol = structConfLoop.find_tag("_struct_conf.id");
        int typeIdCol = structConfLoop.find_tag("_struct_conf.conf_type_id");
        // int beginCompCol = structConfLoop.find_tag("_struct_conf.beg_label_comp_id");
        int beginAsymIdCol = structConfLoop.find_tag("_struct_conf.beg_auth_asym_id");
        int beginLabelSeqCol = structConfLoop.find_tag("_struct_conf.beg_label_seq_id");
        int endLabelSeqCol = structConfLoop.find_tag("_struct_conf.end_label_seq_id");
        for (int row = 0; row < loopLength; row++) {
            const std::string& asymId = structConfLoop.val(row, beginAsymIdCol);
            const int seqBegin = std::stoi(structConfLoop.val(row, beginLabelSeqCol), nullptr);
            const int seqEnd = std::stoi(structConfLoop.val(row, endLabelSeqCol), nullptr);
            const std::string& type = structConfLoop.val(row, typeIdCol);
            SecondaryStruct secondaryStruct = {
                .code = getSecondaryStructureType(type),
                .labelAsymId = asymId,
                .sequenceBegin = seqBegin,
                .sequenceEnd = seqEnd
            };
            if (!secondaryStructures.contains(asymId)) {
                std::vector<SecondaryStruct> structs;
                secondaryStructures[asymId] = structs;
            }
            secondaryStructures[asymId].emplace_back(secondaryStruct);
        }
        return secondaryStructures;
    }

    std::map<string, std::vector<SecondaryStruct>> parseSheetRange(gemmi::cif::Block& block) {
        // Beta strands
        std::map<string, std::vector<SecondaryStruct>> secondaryStructures;
        auto loopPointer = block.find_loop_item("_struct_sheet_range.id");
        if (loopPointer == NULL) {
            return secondaryStructures;
        }
        auto sheetRangeLoop = loopPointer->loop;
        int loopLength = sheetRangeLoop.length();
        // int sheetIdCol = sheetRangeLoop.find_tag("_struct_sheet_range.sheet_id");
        // int idCol = sheetRangeLoop.find_tag("_struct_sheet_range.id");
        // int beginCompCol = sheetRangeLoop.find_tag("_struct_sheet_range.beg_label_comp_id");
        int beginAsymIdCol = sheetRangeLoop.find_tag("_struct_sheet_range.beg_auth_asym_id");
        int beginLabelSeqCol = sheetRangeLoop.find_tag("_struct_sheet_range.beg_label_seq_id");
        int endLabelSeqCol = sheetRangeLoop.find_tag("_struct_sheet_range.end_label_seq_id");
        for (int row = 0; row < loopLength; row++) {
            const std::string& asymId = sheetRangeLoop.val(row, beginAsymIdCol);
            const int seqBegin = std::stoi(sheetRangeLoop.val(row, beginLabelSeqCol), nullptr);
            const int seqEnd = std::stoi(sheetRangeLoop.val(row, endLabelSeqCol));
            SecondaryStruct secondaryStruct = {
                .code = 'E',
                .labelAsymId = asymId,
                .sequenceBegin = seqBegin,
                .sequenceEnd = seqEnd
            };
            if (!secondaryStructures.contains(asymId)) {
                std::vector<SecondaryStruct> structs;
                secondaryStructures[asymId] = structs;
            }
            secondaryStructures[asymId].emplace_back(secondaryStruct);
        }
        return secondaryStructures;
    }

    void updateDsspStrings(std::map<string, string>& dsspStrings, std::map<string, std::vector<SecondaryStruct>>& secondaryStructs) {
        for (const auto& itemPair : secondaryStructs) {
            auto& dssp = dsspStrings[itemPair.first];
            for (const auto& secondaryStruct : itemPair.second) {
                auto length = secondaryStruct.sequenceEnd - secondaryStruct.sequenceBegin + 1;
                auto begin = dssp.begin() + secondaryStruct.sequenceBegin - 1;
                auto end = begin + length;
                for (auto it = begin; it < end; it++) {
                    *it = secondaryStruct.code;
                }
            }
        }
    }
}
