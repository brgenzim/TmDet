#include <string>
#include <vector>
#include <iostream>
#include <format>
#include <pugixml.hpp>
#include <DTOs/Xml/Writer3.hpp>
#include <DTOs/Xml/Constants3.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/TMatrix.hpp>
#include <ValueObjects/Chain.hpp>
#include <Exceptions/SyntaxErrorException.hpp>

namespace Tmdet::DTOs::Xml {

    void Writer3::create() {
        _doc.load_string(_pdbtm_xml.c_str());
        _root = _doc.child(XML_NODE_ROOT);
    }
            
    void Writer3::setTmp(const bool& tmp) const {
        _root.attribute(XML_ATTR_TMP).set_value(tmp?"yes":"no");
    }
    
    void Writer3::setCode(const std::string& code) const {
        _root.attribute(XML_ATTR_ID).set_value(code.c_str());
    }
            
    void Writer3::setCreateDate(const std::string& date) const {
        _root.child(XML_NODE_CREATE_DATE).text() = date.c_str();
    }

    void Writer3::setModifications(const std::vector<Tmdet::ValueObjects::Modification>& mods) {
        pugi::xml_node node;
        for(const auto& m : mods) {
            node = _root.insert_child_after(XML_NODE_MODIFICATION, _root.child(XML_NODE_CREATE_DATE));
            pugi::xml_node date_node = node.append_child(XML_NODE_DATE);
            date_node.append_child(pugi::node_pcdata).set_value(m.date.c_str());
            pugi::xml_node descr_node = node.append_child(XML_NODE_DESCR);
            descr_node.append_child(pugi::node_pcdata).set_value(m.descr.c_str());
        }
    }

    void Writer3::setQvalue(const double& q) const {
        _root.child(XML_NODE_RAWRES).child(XML_NODE_TMRES).text() = std::format("{:.2f}",q).c_str();
    }

    void Writer3::setTmtype(const std::string& ptype) const {
        _root.child(XML_NODE_RAWRES).child(XML_NODE_TMTYPE).text() = ptype.c_str();
    }

    void Writer3::setSpres(const std::string& spres) const {
        _root.child(XML_NODE_RAWRES).child(XML_NODE_SPRES).text() = spres.c_str();
    }

    void Writer3::setPdbkwres(const std::string& pdbkwres) const {
        _root.child(XML_NODE_RAWRES).child(XML_NODE_PDBKWRES).text() = pdbkwres.c_str();
    }

    void Writer3::setTMatrix(const pugi::xml_node& node, Tmdet::ValueObjects::TMatrix& tmatrix) const {
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

    void Writer3::setBioMatrix(Tmdet::ValueObjects::BioMatrix& bioMatrix) {
        if (bioMatrix.matrices.empty() && bioMatrix.deletedChainIds.empty() ) {
            return;
        }
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

    void Writer3::setMembranes(std::vector<Tmdet::ValueObjects::Membrane>& membranes) {
        
        for(auto& m : membranes) {
            pugi::xml_document doc;
            doc.load_string(_membrane_xml.c_str());
            pugi::xml_node node = doc.child(XML_NODE_MEMBRANE);
            node.child(XML_NODE_NORMAL).attribute(XML_ATTR_Z) = std::to_string(m.halfThickness).c_str();
            //setTMatrix(node,m.tmatrix);
            _root.append_copy(node);
        }
    }

    void Writer3::setChains(const std::vector<Tmdet::ValueObjects::Chain>& chains) {
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

    void Writer3::setRegions(pugi::xml_node& pnode, const std::vector<Tmdet::ValueObjects::Region>& regions) const {
        for(const auto& r: regions) {
            pugi::xml_node node = pnode.append_child(XML_NODE_REGION);
            node.append_attribute(XML_ATTR_SEQ_BEG) = std::to_string(r.beg).c_str();
            node.append_attribute(XML_ATTR_PDB_BEG) = r.begi.c_str();
            node.append_attribute(XML_ATTR_SEQ_END) = std::to_string(r.end).c_str();
            node.append_attribute(XML_ATTR_PDB_END) = r.endi.c_str();
            node.append_attribute(XML_ATTR_type) = std::format("{}","%c",r.type.code).c_str();
        }
    }
            
    void Writer3::writeXml(Tmdet::ValueObjects::Protein& protein, const std::string& path) {
        create();
        setTmp(protein.tmp);
        setCode(protein.code);
        setCreateDate(protein.date);
        setModifications(protein.modifications);
        setQvalue(protein.qValue);
        setTmtype(protein.type.name);
        setSpres(protein.spres);
        setPdbkwres(protein.pdbkwres);
        setBioMatrix(protein.bioMatrix);
        setMembranes(protein.membranes);
        setChains(protein.chains);
        write(path);
    }

}
