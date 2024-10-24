#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <gemmi/cifdoc.hpp>
#include <gemmi/model.hpp>
#include <gemmi/modify.hpp>
#include <gemmi/polyheur.hpp>
#include <gemmi/align.hpp>
#include <gemmi/seqalign.hpp>
#include <gemmi/seqtools.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/to_mmcif.hpp>
#include <Config.hpp>
#include <DTOs/Protein.hpp>
#include <System/FilePaths.hpp>
#include <Utils/Alignment.hpp>
#include <Types/Protein.hpp>
#include <Types/Residue.hpp>
#include <Types/SecStruct.hpp>
#include <System/Logger.hpp>
#include <ValueObjects/Protein.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/Residue.hpp>
#include <ValueObjects/Atom.hpp>

namespace Tmdet::DTOs {

    void Protein::writeCif(const Tmdet::ValueObjects::Protein& protein, const std::string& path) {
        std::ofstream outCif(path);
        gemmi::cif::WriteOptions options(gemmi::cif::Style::Pdbx);
        gemmi::cif::Document document = make_mmcif_document(protein.gemmi);

        // correction of _chem_comp
        logger.debug("Updating chem_comp types before writing document into '{}'", path);
        auto& newChemLoop = document.blocks[0].init_mmcif_loop("_chem_comp.", { "id", "type" });
        auto oldBlock = protein.document.blocks[0];
        for (auto chemComp : oldBlock.find("_chem_comp.", { "id", "type" })) {
            logger.debug("chem_comp: {} {}", chemComp[0], chemComp[1]);
            newChemLoop.add_row({ chemComp[0], chemComp[1] });
        }

        gemmi::cif::write_cif_to_stream(outCif, document, options);
    }

    Tmdet::ValueObjects::Protein Protein::get(const std::string& inputPath) {
        logger.debug("Processing Protein::get()");
        Tmdet::ValueObjects::Protein protein;
        protein.getStructure(inputPath);
        protein.code = protein.gemmi.name;
        remove_hydrogens(protein.gemmi.models[0]);
        remove_ligands_and_waters(protein.gemmi.models[0]);
        remove_alternative_conformations(protein.gemmi.models[0]);
        // keep only the first model
        protein.gemmi.models.resize(1);

        // Fill residue gaps in chains where there is usable entity sequence data
        Tmdet::Utils::Alignment::alignResidues(protein);
        int chainIdx = 0;
        for(auto& chain: protein.gemmi.models[0].chains) {
            Tmdet::ValueObjects::Chain chainVO;
            chainVO.addStructure(chain);
            chainVO.idx = chainIdx++;
            auto poly = chain.get_polymer();
            chainVO.labId = (poly?poly.subchain_id():chainVO.id);
            logger.debug("Loading chain auth_asym_id: {} label_asym_id: {}",chainVO.id,chainVO.labId);
            int residueIdx = 0;
            for( auto& residue: chain.residues) {
                auto residueVO = Tmdet::ValueObjects::Residue(residue,chainVO);
                residueVO.chainIdx = chainVO.idx;
                residueVO.idx = residueIdx++;
                residueVO.surface = 0.0;
                chainVO.seq += residueVO.type.a1;
                int atomIdx = 0;
                for( auto& atom: residue.atoms) {
                    auto atomVO = Tmdet::ValueObjects::Atom(atom,residueVO,chainVO);
                    atomVO.chainIdx = chainVO.idx;
                    atomVO.residueIdx = residueVO.idx;
                    atomVO.idx = atomIdx++;
                    atomVO.surface = 0.0;
                    residueVO.atoms.emplace_back(atomVO);
                }
                chainVO.residues.emplace_back(residueVO);
            }
            chainVO.length = residueIdx;
            protein.chains.emplace_back(chainVO);
        }

        protein.neighbors = gemmi::NeighborSearch(protein.gemmi.models[0], protein.gemmi.cell, 9);
        protein.neighbors.populate();
        logger.debug(" Processed Protein::get()");
        return protein;
    }

    void Protein::unselectPolymers(Tmdet::ValueObjects::Protein& protein) {
        // Load unselectable names
        auto filterPath = environment.get("TMDET_POLYMER_FILTER_FILE", DEFAULT_TMDET_POLYMER_FILTER_FILE);
        if ( !Tmdet::System::FilePaths::fileExists(filterPath) ) {
            logger.warn("polymer filter file not found: {}", filterPath);
            return;
        }

        std::ifstream filters;
        filters.open(filterPath);
        std::set<std::string> nameSet;
        for (std::string line; std::getline(filters, line); ) {
            if (line.length() == 0) {
                continue;
            }
            nameSet.insert(line);
        }
        filters.close();

        // to_upper implementation
        auto toUpper = [](std::string subject) {
            std::transform(subject.begin(), subject.end(), subject.begin(),
                [](unsigned char c){ return std::toupper(c); }
            );
            return subject;
        };

        for (auto& chain : protein.chains) {
            if (!protein.polymerNames.contains(chain.entityId)) {
                continue;
            }
            auto name = protein.polymerNames[chain.entityId];
            for (auto& filter : nameSet) {
                if (toUpper(name).find(toUpper(filter)) != name.npos) {
                    chain.selected = false;
                    break;
                }
            }
        }
    }

    void Protein::print(std::ostream& os, const Tmdet::ValueObjects::Protein& protein) {
        for(const auto& chainVO: protein.chains) {
            printChain(os, chainVO);
        }
    }

    void Protein::printChain(std::ostream& os, const Tmdet::ValueObjects::Chain& chainVO) {
        for(const auto& residueVO: chainVO.residues) {
            os << "CHAIN " << chainVO.gemmi.name << " " << chainVO.length << std::endl;
            printResidue(os, residueVO);
        }
    }

    void Protein::printResidue(std::ostream& os, const Tmdet::ValueObjects::Residue& residueVO) {
        os << "\tRESIDUE " << residueVO.idx << ":" << residueVO.resn() << "(" << residueVO.gemmi.name << ") ";
        os << residueVO.surface << " " << residueVO.ss.code << std::endl;
        if ( residueVO.temp.contains("fragment") ) {
            os << "\t\tTEMP: fragment: " << any_cast<int>(residueVO.temp.at("fragment")) << std::endl;
        }
        for(const auto& atomVO: residueVO.atoms) {
            printAtom(os, atomVO);
        }
    }

    void Protein::printAtom(std::ostream& os, const Tmdet::ValueObjects::Atom& atomVO) {
        os << "\t\tATOM " << atomVO.idx << ": " << atomVO.gemmi.name << " ";
        os << atomVO.gemmi.pos.x << " " << atomVO.gemmi.pos.y << " " << atomVO.gemmi.pos.z << " ";
        Tmdet::Types::Residue residueType = Tmdet::Types::ResidueType::getResidue(atomVO.parentResidue.gemmi.name);
        double vdw = 0.0;
        if (residueType.atoms.contains(atomVO.gemmi.name)) {
            vdw = residueType.atoms.at(atomVO.gemmi.name).atom.vdw;
        } else {
            vdw = Types::AtomType::DEFAULT_VDW;
        }
        os << atomVO.surface << " " << vdw;
        if (atomVO.temp.contains("outside")) {
            os << " out: " << any_cast<double>(atomVO.temp.at("outside"));
        }
        os << std::endl;
    }

    std::vector<std::string> Protein::getChainSequence(
        const Tmdet::ValueObjects::Protein& protein, const gemmi::Chain& chain) {

        std::vector<std::string> sequence;
        auto entityId = chain.residues[0].entity_id;
        for (const auto& entity : protein.gemmi.entities) {
            if (entity.entity_type == gemmi::EntityType::Polymer
                && entity.name == entityId) {

                sequence = entity.full_sequence;
                break;
            }
        }
        return sequence;
    }

    void Protein::transform(Tmdet::ValueObjects::Protein& protein) {
        for (auto& chain: protein.chains) {
            for (auto& residue: chain.residues) {
                for (auto& atom: residue.atoms) {
                    double x = protein.tmatrix.trans.x
                                + atom.gemmi.pos.x * protein.tmatrix.rot[0][0]
                                + atom.gemmi.pos.y * protein.tmatrix.rot[0][1]
                                + atom.gemmi.pos.z * protein.tmatrix.rot[0][2];
			        double y = protein.tmatrix.trans.y
                                + atom.gemmi.pos.x * protein.tmatrix.rot[1][0]
                                + atom.gemmi.pos.y * protein.tmatrix.rot[1][1]
                                + atom.gemmi.pos.z * protein.tmatrix.rot[1][2];
			        double z = protein.tmatrix.trans.z
                                + atom.gemmi.pos.x * protein.tmatrix.rot[2][0]
                                + atom.gemmi.pos.y * protein.tmatrix.rot[2][1]
                                + atom.gemmi.pos.z * protein.tmatrix.rot[2][2];
                    atom.gemmi.pos.x = x;
                    atom.gemmi.pos.y = y;
                    atom.gemmi.pos.z = z;
                }
            }
        }
    }
}
