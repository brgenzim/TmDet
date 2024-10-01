#include <string>
#include <vector>
#include <iostream>
#include <format>
#include <pugixml.hpp>
#include <DTOs/Xml.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/TMatrix.hpp>
#include <ValueObjects/Chain.hpp>
#include <Exceptions/SyntaxErrorException.hpp>

namespace Tmdet::DTOs {
    
    void Xml::read(const std::string& path) {
        if (pugi::xml_parse_result result = _doc.load_file(path.c_str()); !result) {
            throw Tmdet::Exceptions::SyntaxErrorException(path,(int)result.offset,result.description());
        }
        _root = _doc.child(XML_NODE_ROOT);
    }

    void Xml::write(const std::string& path) const {
        _doc.save_file(path.c_str(),"  ");
    }

    void Xml::create() {
        _doc.load_string(_pdbtm_xml.c_str());
        _root = _doc.child(XML_NODE_ROOT);
    }

    bool Xml::getTmp() const {
        return _root.attribute(XML_ATTR_TMP).as_bool();
    }

    void Xml::setTmp(const bool& tmp) const {
        _root.attribute(XML_ATTR_TMP).set_value(tmp?"yes":"no");
    }

    std::string Xml::getCode() const {
        return _root.attribute(XML_ATTR_ID).value();
    }

    void Xml::setCode(const std::string& code) const {
        _root.attribute(XML_ATTR_ID).set_value(code.c_str());
    }

    std::string Xml::getCreateDate() const {
        return _root.child(XML_NODE_CREATE_DATE).text().get();
    }

    void Xml::setCreateDate(const std::string& date) const {
        _root.child(XML_NODE_CREATE_DATE).text() = date.c_str();
    }

    std::vector<Tmdet::ValueObjects::Modification> Xml::getModifications() const {
        std::vector<Tmdet::ValueObjects::Modification> mods;
        for (pugi::xml_node mod = _root.child(XML_NODE_MODIFICATION); mod; mod = mod.next_sibling(XML_NODE_MODIFICATION)) {
            mods.emplace_back(mod.child(XML_NODE_DATE).text().get(),mod.child(XML_NODE_DESCR).text().get());
        }
        return mods;
    }

    void Xml::setModifications(const std::vector<Tmdet::ValueObjects::Modification>& mods) {
        pugi::xml_node node;
        for(const auto& m : mods) {
            node = _root.insert_child_after(XML_NODE_MODIFICATION, _root.child(XML_NODE_CREATE_DATE));
            pugi::xml_node date_node = node.append_child(XML_NODE_DATE);
            date_node.append_child(pugi::node_pcdata).set_value(m.date.c_str());
            pugi::xml_node descr_node = node.append_child(XML_NODE_DESCR);
            descr_node.append_child(pugi::node_pcdata).set_value(m.descr.c_str());
        }
    }

    double Xml::getQvalue() const {
        return _root.child(XML_NODE_RAWRES).child(XML_NODE_TMRES).text().as_double();
    }

    void Xml::setQvalue(const double& q) const {
        _root.child(XML_NODE_RAWRES).child(XML_NODE_TMRES).text() = std::format("{}","%.2f",q).c_str();
    }

    std::string Xml::getTmtype() const {
        return _root.child(XML_NODE_RAWRES).child(XML_NODE_TMTYPE).text().get();
    }

    void Xml::setTmtype(const std::string& ptype) const {
        _root.child(XML_NODE_RAWRES).child(XML_NODE_TMTYPE).text() = ptype.c_str();
    }

    std::string Xml::getSpres() const {
        return _root.child(XML_NODE_RAWRES).child(XML_NODE_SPRES).text().get();
    }

    void Xml::setSpres(const std::string& spres) const {
        _root.child(XML_NODE_RAWRES).child(XML_NODE_SPRES).text() = spres.c_str();
    }

    std::string Xml::getPdbkwres() const {
        return _root.child(XML_NODE_RAWRES).child(XML_NODE_PDBKWRES).text().get();
    }

    void Xml::setPdbkwres(const std::string& pdbkwres) const {
        _root.child(XML_NODE_RAWRES).child(XML_NODE_PDBKWRES).text() = pdbkwres.c_str();
    }

    Tmdet::ValueObjects::TMatrix Xml::getTMatrix(const pugi::xml_node& node) const {
        Tmdet::ValueObjects::TMatrix tmatrix;
        tmatrix.rot[0][0] = node.child(XML_NODE_ROWX).attribute(XML_ATTR_X).as_double();
        tmatrix.rot[0][1] = node.child(XML_NODE_ROWX).attribute(XML_ATTR_Y).as_double();
        tmatrix.rot[0][2] = node.child(XML_NODE_ROWX).attribute(XML_ATTR_Z).as_double();
        tmatrix.rot[1][0] = node.child(XML_NODE_ROWY).attribute(XML_ATTR_X).as_double();
        tmatrix.rot[1][1] = node.child(XML_NODE_ROWY).attribute(XML_ATTR_Y).as_double();
        tmatrix.rot[1][2] = node.child(XML_NODE_ROWY).attribute(XML_ATTR_Z).as_double();
        tmatrix.rot[2][0] = node.child(XML_NODE_ROWZ).attribute(XML_ATTR_X).as_double();
        tmatrix.rot[2][1] = node.child(XML_NODE_ROWZ).attribute(XML_ATTR_Y).as_double();
        tmatrix.rot[2][2] = node.child(XML_NODE_ROWZ).attribute(XML_ATTR_Z).as_double();
        tmatrix.trans.x = node.child(XML_NODE_ROWX).attribute(XML_ATTR_T).as_double();
        tmatrix.trans.y = node.child(XML_NODE_ROWY).attribute(XML_ATTR_T).as_double();
        tmatrix.trans.z = node.child(XML_NODE_ROWZ).attribute(XML_ATTR_T).as_double();
        return tmatrix;
    }

    void Xml::setTMatrix(const pugi::xml_node& node, Tmdet::ValueObjects::TMatrix& tmatrix) const {
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWX).attribute(XML_ATTR_X) = std::to_string(tmatrix.rot[0][0]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWX).attribute(XML_ATTR_Y) = std::to_string(tmatrix.rot[0][1]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWX).attribute(XML_ATTR_Z) = std::to_string(tmatrix.rot[0][2]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWY).attribute(XML_ATTR_X) = std::to_string(tmatrix.rot[1][0]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWY).attribute(XML_ATTR_Y) = std::to_string(tmatrix.rot[1][1]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWY).attribute(XML_ATTR_Z) = std::to_string(tmatrix.rot[1][2]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWZ).attribute(XML_ATTR_X) = std::to_string(tmatrix.rot[2][0]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWZ).attribute(XML_ATTR_Y) = std::to_string(tmatrix.rot[2][1]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWZ).attribute(XML_ATTR_Z) = std::to_string(tmatrix.rot[2][2]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWX).attribute(XML_ATTR_T) = std::to_string(tmatrix.trans.x).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWY).attribute(XML_ATTR_T) = std::to_string(tmatrix.trans.y).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWZ).attribute(XML_ATTR_T) = std::to_string(tmatrix.trans.z).c_str();
    }

     Tmdet::ValueObjects::BioMatrix Xml::getBioMatrix() const {
        Tmdet::ValueObjects::BioMatrix bioMatrix;
        pugi::xml_node node = _root.child(XML_NODE_BIOMATRIX);
        for (pugi::xml_node matrix = node.child(XML_NODE_MATRIX); matrix; matrix = matrix.next_sibling(XML_NODE_MATRIX)) {
            pugi::xml_node tnode = matrix.child(XML_NODE_TMATRIX);
            bioMatrix.matrices.emplace_back(
                matrix.attribute(XML_ATTR_ID).as_int(),
                matrix.child(XML_NODE_APPLY_TO_CHAIN).attribute(XML_ATTR_CHAINID).as_string(),
                matrix.child(XML_NODE_APPLY_TO_CHAIN).attribute(XML_ATTR_NEW_CHAINID).as_string(),
                getTMatrix(tnode)
            );
        }
        return bioMatrix;
     }

    void Xml::setBioMatrix(Tmdet::ValueObjects::BioMatrix& bioMatrix) {
        pugi::xml_node biom_node = _root.insert_child_after(XML_NODE_BIOMATRIX, _root.child(XML_NODE_RAWRES));
        for(auto& m: bioMatrix.matrices) {
            pugi::xml_document doc;
            doc.load_string(_biomatrix_xml.c_str());
            pugi::xml_node node = doc.child(XML_NODE_MATRIX);
            node.attribute(XML_ATTR_ID) = std::to_string(m.id).c_str();
            node.child(XML_NODE_APPLY_TO_CHAIN).attribute(XML_ATTR_CHAINID) = m.sourceChainId.c_str();
            node.child(XML_NODE_APPLY_TO_CHAIN).attribute(XML_ATTR_NEW_CHAINID) = m.newChainId.c_str();
            setTMatrix(node, m.tmatrix);
            biom_node.append_copy(node);
        }

    }

    std::vector<Tmdet::ValueObjects::Membrane> Xml::getMembranes() const {
        std::vector<Tmdet::ValueObjects::Membrane> membranes;
        for (pugi::xml_node m_node = _root.child(XML_NODE_MEMBRANE); m_node; m_node = m_node.next_sibling(XML_NODE_MEMBRANE)) {
            pugi::xml_node tnode = m_node.child(XML_NODE_TMATRIX);
            membranes.emplace_back(getTMatrix(tnode),
                gemmi::Vec3(0,0,0), //todo get origo
                gemmi::Vec3(0,0,1), //todo get normal
                m_node.child(XML_NODE_NORMAL).attribute(XML_ATTR_Z).as_double(),
                0.0,
                0.0,
                Tmdet::Types::Membranes.at("Plain")
            );
        }
        return membranes;
    }

    void Xml::setMembranes(std::vector<Tmdet::ValueObjects::Membrane>& membranes) {
        
        for(auto& m : membranes) {
            pugi::xml_document doc;
            doc.load_string(_membrane_xml.c_str());
            pugi::xml_node node = doc.child(XML_NODE_MEMBRANE);
            node.child(XML_NODE_NORMAL).attribute(XML_ATTR_Z) = std::to_string(m.h).c_str();
            setTMatrix(node,m.tmatrix);
            _root.append_copy(node);
        }
    }

    void Xml::getChains(std::vector<Tmdet::ValueObjects::Chain>& chains) {
        for (pugi::xml_node c_node = _root.child(XML_NODE_CHAIN); c_node; c_node = c_node.next_sibling(XML_NODE_CHAIN)) {
            std::string type = c_node.attribute(XML_ATTR_TYPE).as_string();
            bool found = false;
            for( auto& c: chains) {
                if (c.id == c_node.attribute(XML_ATTR_CHAINID).as_string()) {
                    c.selected = true;
                    c.numtm = c_node.attribute(XML_ATTR_NUM_TM).as_int();
                    c.seq = c_node.child(XML_NODE_SEQ).text().get();
                    c.regions = getRegions(c_node);
                    c.type = Tmdet::Types::Chains.at(type);
                    found = true;
                    continue;
                }
            }
            if (!found) {
                std::cerr << "Could not find gemmi chain for " << c_node.attribute(XML_ATTR_CHAINID).as_string() << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    void Xml::setChains(const std::vector<Tmdet::ValueObjects::Chain>& chains) {
        for(const auto& c: chains) {
            pugi::xml_document doc;
            doc.load_string(_chain_xml.c_str());
            pugi::xml_node node = doc.child(XML_NODE_CHAIN);
            node.attribute(XML_ATTR_CHAINID) = c.id.c_str();
            node.attribute(XML_ATTR_NUM_TM) = std::to_string(c.numtm).c_str();
            node.attribute(XML_ATTR_TYPE) = c.type.name.c_str();
            node.child(XML_NODE_SEQ).text() = c.seq.c_str();
            setRegions(node,c.regions);
            _root.append_copy(node);
        }
    }

    std::vector<Tmdet::ValueObjects::Region> Xml::getRegions(const pugi::xml_node& cnode) const {
        std::vector<Tmdet::ValueObjects::Region> regions;
        for (pugi::xml_node r_node = cnode.child(XML_NODE_REGION); r_node; r_node = r_node.next_sibling(XML_NODE_REGION)) {
            char type = r_node.attribute(XML_ATTR_type).value()[0];
            regions.emplace_back(r_node.attribute(XML_ATTR_SEQ_BEG).as_int(),
                r_node.attribute(XML_ATTR_PDB_BEG).as_string(),
                r_node.attribute(XML_ATTR_SEQ_END).as_int(),
                r_node.attribute(XML_ATTR_PDB_END).as_string(),
                0,
                0,
                Tmdet::Types::Regions.at(type));
        }
        return regions;
    }

    void Xml::setRegions(pugi::xml_node& pnode, const std::vector<Tmdet::ValueObjects::Region>& regions) const {
        for(const auto& r: regions) {
            pugi::xml_node node = pnode.append_child(XML_NODE_REGION);
            node.append_attribute(XML_ATTR_SEQ_BEG) = std::to_string(r.beg).c_str();
            node.append_attribute(XML_ATTR_PDB_BEG) = r.begi.c_str();
            node.append_attribute(XML_ATTR_SEQ_END) = std::to_string(r.end).c_str();
            node.append_attribute(XML_ATTR_PDB_END) = r.endi.c_str();
            node.append_attribute(XML_ATTR_type) = std::format("{}","%c",r.type.code).c_str();
        }
    }

}