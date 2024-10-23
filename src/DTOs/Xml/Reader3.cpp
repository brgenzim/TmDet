#include <string>
#include <vector>
#include <iostream>
#include <format>
#include <pugixml.hpp>
#include <DTOs/Xml/Constants3.hpp>
#include <DTOs/Xml/Reader3.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/TMatrix.hpp>
#include <ValueObjects/Chain.hpp>
#include <Exceptions/SyntaxErrorException.hpp>

namespace Tmdet::DTOs::Xml {

    void Reader3::read(const std::string& path) {
        if (pugi::xml_parse_result result = _doc.load_file(path.c_str()); !result) {
            throw Tmdet::Exceptions::SyntaxErrorException(path,(int)result.offset,result.description());
        }
        _root = _doc.child(XML3_NODE_ROOT);
    }

    void Reader3::readXml(Tmdet::ValueObjects::Protein& protein, const std::string& path) {
        read(path);
        protein.tmp = getTmp();
        protein.code = getCode();
        protein.date = getCreateDate();
        protein.modifications = getModifications();
        protein.qValue = getQvalue();
        protein.type = Tmdet::Types::Proteins.at(getTmtype());
        protein.spres = getSpres();
        protein.pdbkwres = getPdbkwres();
        protein.bioMatrix = getBioMatrix();
        protein.membranes = getMembranes();
        getChains(protein.chains);
    }

    bool Reader3::getTmp() const {
        return _root.attribute(XML3_ATTR_TMP).as_bool();
    }

    std::string Reader3::getCode() const {
        return _root.attribute(XML3_ATTR_ID).value();
    }

    std::string Reader3::getCreateDate() const {
        return _root.child(XML3_NODE_CREATE_DATE).text().get();
    }

    std::vector<Tmdet::ValueObjects::Modification> Reader3::getModifications() const {
        std::vector<Tmdet::ValueObjects::Modification> mods;
        for (pugi::xml_node mod = _root.child(XML3_NODE_MODIFICATION); mod; mod = mod.next_sibling(XML3_NODE_MODIFICATION)) {
            mods.emplace_back(mod.child(XML3_NODE_DATE).text().get(),mod.child(XML3_NODE_DESCR).text().get());
        }
        return mods;
    }

    double Reader3::getQvalue() const {
        return _root.child(XML3_NODE_RAWRES).child(XML3_NODE_TMRES).text().as_double();
    }

    std::string Reader3::getTmtype() const {
        return _root.child(XML3_NODE_RAWRES).child(XML3_NODE_TMTYPE).text().get();
    }

    std::string Reader3::getSpres() const {
        return _root.child(XML3_NODE_RAWRES).child(XML3_NODE_SPRES)?
                    _root.child(XML3_NODE_RAWRES).child(XML3_NODE_SPRES).text().get():(std::string)"";
    }

    std::string Reader3::getPdbkwres() const {
        return _root.child(XML3_NODE_RAWRES).child(XML3_NODE_PDBKWRES)?
                    _root.child(XML3_NODE_RAWRES).child(XML3_NODE_PDBKWRES).text().get():(std::string)"";
    }

    Tmdet::ValueObjects::TMatrix Reader3::getTMatrix(const pugi::xml_node& node) const {
        Tmdet::ValueObjects::TMatrix tmatrix;
        tmatrix.rot[0][0] = node.child(XML3_NODE_ROWX).attribute(XML3_ATTR_X).as_double();
        tmatrix.rot[0][1] = node.child(XML3_NODE_ROWX).attribute(XML3_ATTR_Y).as_double();
        tmatrix.rot[0][2] = node.child(XML3_NODE_ROWX).attribute(XML3_ATTR_Z).as_double();
        tmatrix.rot[1][0] = node.child(XML3_NODE_ROWY).attribute(XML3_ATTR_X).as_double();
        tmatrix.rot[1][1] = node.child(XML3_NODE_ROWY).attribute(XML3_ATTR_Y).as_double();
        tmatrix.rot[1][2] = node.child(XML3_NODE_ROWY).attribute(XML3_ATTR_Z).as_double();
        tmatrix.rot[2][0] = node.child(XML3_NODE_ROWZ).attribute(XML3_ATTR_X).as_double();
        tmatrix.rot[2][1] = node.child(XML3_NODE_ROWZ).attribute(XML3_ATTR_Y).as_double();
        tmatrix.rot[2][2] = node.child(XML3_NODE_ROWZ).attribute(XML3_ATTR_Z).as_double();
        tmatrix.trans.x = node.child(XML3_NODE_ROWX).attribute(XML3_ATTR_T).as_double();
        tmatrix.trans.y = node.child(XML3_NODE_ROWY).attribute(XML3_ATTR_T).as_double();
        tmatrix.trans.z = node.child(XML3_NODE_ROWZ).attribute(XML3_ATTR_T).as_double();
        return tmatrix;
    }

    Tmdet::ValueObjects::BioMatrix Reader3::getBioMatrix() const {
        Tmdet::ValueObjects::BioMatrix bioMatrix;
        pugi::xml_node node = _root.child(XML3_NODE_BIOMATRIX);
        for (pugi::xml_node matrix = node.child(XML3_NODE_MATRIX); matrix; matrix = matrix.next_sibling(XML3_NODE_MATRIX)) {
            pugi::xml_node tnode = matrix.child(XML3_NODE_TMATRIX);
            bioMatrix.matrices.emplace_back(
                matrix.attribute(XML3_ATTR_ID).as_int(),
                matrix.child(XML3_NODE_APPLY_TO_CHAIN).attribute(XML3_ATTR_CHAINID).as_string(),
                matrix.child(XML3_NODE_APPLY_TO_CHAIN).attribute(XML3_ATTR_NEW_CHAINID).as_string(),
                getTMatrix(tnode)
            );
        }
        return bioMatrix;
    }

    std::vector<Tmdet::ValueObjects::Membrane> Reader3::getMembranes() const {
        std::vector<Tmdet::ValueObjects::Membrane> membranes;
        for (pugi::xml_node m_node = _root.child(XML3_NODE_MEMBRANE); m_node; m_node = m_node.next_sibling(XML3_NODE_MEMBRANE)) {
            pugi::xml_node tnode = m_node.child(XML3_NODE_TMATRIX);
            membranes.emplace_back(//getTMatrix(tnode),
                0.0, //todo get origo
                m_node.child(XML3_NODE_NORMAL).attribute(XML3_ATTR_Z).as_double(),
                0.0,
                10.0,
                Tmdet::Types::Membranes.at("Plain")
            );
        }
        return membranes;
    }

    void Reader3::getChains(std::vector<Tmdet::ValueObjects::Chain>& chains) {
        for (pugi::xml_node c_node = _root.child(XML3_NODE_CHAIN); c_node; c_node = c_node.next_sibling(XML3_NODE_CHAIN)) {
            std::string type = c_node.attribute(XML3_ATTR_TYPE).as_string();
            bool found = false;
            for( auto& c: chains) {
                if (c.id == c_node.attribute(XML3_ATTR_CHAINID).as_string()) {
                    c.selected = true;
                    c.numtm = c_node.attribute(XML3_ATTR_NUM_TM).as_int();
                    c.seq = c_node.child(XML3_NODE_SEQ).text().get();
                    c.regions = getRegions(c_node);
                    c.type = Tmdet::Types::Chains.at(type);
                    found = true;
                    continue;
                }
            }
            if (!found) {
                Tmdet::ValueObjects::Chain c;
                c.id = c_node.attribute(XML3_ATTR_CHAINID).as_string();
                c.selected = true;
                c.numtm = c_node.attribute(XML3_ATTR_NUM_TM).as_int();
                c.seq = c_node.child(XML3_NODE_SEQ).text().get();
                c.regions = getRegions(c_node);
                c.type = Tmdet::Types::Chains.at(type);
                chains.emplace_back(c);
            }
        }
    }

    std::vector<Tmdet::ValueObjects::Region> Reader3::getRegions(const pugi::xml_node& cnode) const {
        std::vector<Tmdet::ValueObjects::Region> regions;
        for (pugi::xml_node r_node = cnode.child(XML3_NODE_REGION); r_node; r_node = r_node.next_sibling(XML3_NODE_REGION)) {
            char type = r_node.attribute(XML3_ATTR_type).value()[0];
            regions.emplace_back(r_node.attribute(XML3_ATTR_SEQ_BEG).as_int(),
                r_node.attribute(XML3_ATTR_PDB_BEG).as_int(),
                (r_node.attribute(XML3_ATTR_PDB_BEGI)?r_node.attribute(XML3_ATTR_PDB_BEGI).as_string()[0]:' '),
                r_node.attribute(XML3_ATTR_PDB_BEG).as_int(),
                r_node.attribute(XML3_ATTR_SEQ_END).as_int(),
                r_node.attribute(XML3_ATTR_PDB_END).as_int(),
                (r_node.attribute(XML3_ATTR_PDB_ENDI)?r_node.attribute(XML3_ATTR_PDB_ENDI).as_string()[0]:' '),
                r_node.attribute(XML3_ATTR_PDB_END).as_int(),
                Tmdet::Types::Regions.at(type));
        }
        return regions;
    }

}
