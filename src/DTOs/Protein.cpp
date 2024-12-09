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

    void Protein::writeCif(const Tmdet::VOs::Protein& protein, const std::string& path) {
        gemmi::cif::Document document = make_mmcif_document(protein.gemmi);

        // correction of _chem_comp
        DEBUG_LOG("Updating chem_comp types before writing document into '{}'", path);
        auto& newChemLoop = document.blocks[0].init_mmcif_loop("_chem_comp.", { "id", "type" });
        auto oldBlock = protein.document.blocks[0];
        for (auto chemComp : oldBlock.find("_chem_comp.", { "id", "type" })) {
            DEBUG_LOG("chem_comp: {} {}", chemComp[0], chemComp[1]);
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
