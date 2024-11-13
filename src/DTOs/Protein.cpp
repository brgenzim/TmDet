#include <sstream>
#include <gemmi/modify.hpp>
#include <gemmi/polyheur.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/to_mmcif.hpp>
#include <DTOs/Protein.hpp>
#include <Helpers/Gzip.hpp>
#include <Helpers/String.hpp>
#include <System/FilePaths.hpp>
#include <Utils/Alignment.hpp>
#include <System/Logger.hpp>
#include <ValueObjects/Protein.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/Residue.hpp>
#include <ValueObjects/Atom.hpp>

namespace Tmdet::DTOs {

    void Protein::writeCif(const Tmdet::ValueObjects::Protein& protein, const std::string& path) {
        gemmi::cif::Document document = make_mmcif_document(protein.gemmi);

        // correction of _chem_comp
        logger.debug("Updating chem_comp types before writing document into '{}'", path);
        auto& newChemLoop = document.blocks[0].init_mmcif_loop("_chem_comp.", { "id", "type" });
        auto oldBlock = protein.document.blocks[0];
        for (auto chemComp : oldBlock.find("_chem_comp.", { "id", "type" })) {
            logger.debug("chem_comp: {} {}", chemComp[0], chemComp[1]);
            newChemLoop.add_row({ chemComp[0], chemComp[1] });
        }

        std::stringstream sstream;
        gemmi::cif::WriteOptions options(gemmi::cif::Style::Pdbx);
        gemmi::cif::write_cif_to_stream(sstream, document, options);

        if (path.ends_with(".gz")) {
            Tmdet::Helpers::Gzip::writeFile(path, sstream.str());
        } else {
            std::ofstream outCif(path);
            outCif << sstream.str();
        }
    }

    Tmdet::ValueObjects::Protein Protein::get(const std::string& inputPath) {
        logger.debug("Processing Protein::get()");
        Tmdet::ValueObjects::Protein protein;
        protein.getStructure(inputPath);
        protein.code = protein.gemmi.name;
        Tmdet::Helpers::String::toLower(protein.code);
        remove_hydrogens(protein.gemmi.models[0]);
        remove_ligands_and_waters(protein.gemmi.models[0]);
        remove_alternative_conformations(protein.gemmi.models[0]);
        // keep only the first model
        protein.gemmi.models.resize(1);

        // Fill residue gaps in chains where there is usable entity sequence data
        Tmdet::Utils::Alignment::alignResidues(protein);
        int chainIdx = 0;
        for(auto& chain: protein.gemmi.models[0].chains) {
            if (auto poly = chain.get_polymer(); poly) {
                Tmdet::ValueObjects::Chain chainVO;
                chainVO.addStructure(chain);
                chainVO.idx = chainIdx++;
                chainVO.labId = poly.subchain_id();
                logger.debug("Loading chain auth_asym_id: {} label_asym_id: {}",chainVO.id,chainVO.labId);
                int residueIdx = 0;
                for( auto& residue: chain.residues) {
                    auto residueVO = Tmdet::ValueObjects::Residue(residue);
                    residueVO.chainIdx = chainVO.idx;
                    residueVO.idx = residueIdx++;
                    residueVO.surface = 0.0;
                    chainVO.seq += residueVO.type.a1;
                    int atomIdx = 0;
                    for( auto& atom: residue.atoms) {
                        auto atomVO = Tmdet::ValueObjects::Atom(atom);
                        atomVO.chainIdx = chainVO.idx;
                        atomVO.residueIdx = residueVO.idx;
                        atomVO.idx = atomIdx++;
                        atomVO.surface = 0.0;
                        residueVO.atoms.emplace_back(atomVO);
                    }
                    residueVO.setNumberOfAtoms();
                    chainVO.residues.push_back(residueVO);
                }
                chainVO.length = residueIdx;
                protein.chains.push_back(chainVO);
            }
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
                    DEBUG_LOG("Unselecting Ab chain: {}",chain.id);
                    chain.selected = false;
                    break;
                }
            }
        }
    }

    void Protein::unselectChains(const std::string& chainList, Tmdet::ValueObjects::Protein& protein) {
        for(auto chain: Tmdet::Helpers::String::explode(",",chainList)) {
            if (int chainIdx = protein.searchChainById(chain); chainIdx != -1) {
                protein.chains[chainIdx].selected = false;
                DEBUG_LOG("Unselecting chain: {}",chain);
            }
        }
    }

    void Protein::print(std::ostream& os, const Tmdet::ValueObjects::Protein& protein) {
        for(const auto& chainVO: protein.chains) {
            printChain(os, chainVO);
        }
    }

    void Protein::printChain(std::ostream& os, const Tmdet::ValueObjects::Chain& chainVO) {
        os << "CHAIN " << chainVO.gemmi.name << " " << chainVO.length << std::endl;
        for(const auto& residueVO: chainVO.residues) {
            printResidue(os, residueVO);
        }
    }

    void Protein::printResidue(std::ostream& os, const Tmdet::ValueObjects::Residue& residueVO) {
        os << "\tRESIDUE " << residueVO.idx << ":" << residueVO.resn() << "(" << residueVO.gemmi.name << ") ";
        os << residueVO.surface << " ss:" << residueVO.ss.code << std::endl;
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
            transformChain(chain, protein.tmatrix);
        }
        for (auto& vector: protein.vectors) {
            transformSecStrVec(vector, protein.tmatrix);
        }
    }

    void Protein::transformChain(Tmdet::ValueObjects::Chain& chain, Tmdet::ValueObjects::TMatrix& tmatrix) {
        for (auto& residue: chain.residues) {
            for (auto& atom: residue.atoms) {
                transform(atom.gemmi.pos, tmatrix);
            }
        }
    }

    void Protein::transformSecStrVec(Tmdet::ValueObjects::SecStrVec& vector, Tmdet::ValueObjects::TMatrix& tmatrix) {
        transform(vector.begin, tmatrix);
        transform(vector.end, tmatrix);
    }

    void Protein::transform(gemmi::Vec3& vec, Tmdet::ValueObjects::TMatrix& tmatrix) {
        vec.x += tmatrix.trans.x;
        vec.y += tmatrix.trans.y;
        vec.z += tmatrix.trans.z;
        double x = vec.x * tmatrix.rot[0][0]
                    + vec.y * tmatrix.rot[0][1]
                    + vec.z * tmatrix.rot[0][2];
        double y = vec.x * tmatrix.rot[1][0]
                    + vec.y * tmatrix.rot[1][1]
                    + vec.z * tmatrix.rot[1][2];
        double z = vec.x * tmatrix.rot[2][0]
                    + vec.y * tmatrix.rot[2][1]
                    + vec.z * tmatrix.rot[2][2];
        vec.x = x;
        vec.y = y;
        vec.z = z;
    }
}
