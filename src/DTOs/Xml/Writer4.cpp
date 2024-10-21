#include <string>
#include <vector>
#include <iostream>
#include <format>
#include <pugixml.hpp>
#include <DTOs/Xml/Writer4.hpp>
#include <DTOs/Xml/Constants4.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/TMatrix.hpp>
#include <ValueObjects/Chain.hpp>
#include <Exceptions/SyntaxErrorException.hpp>

namespace Tmdet::DTOs::Xml {

    void Writer4::create() {
        _doc.load_string(_pdbtm_xml.c_str());
        _root = _doc.child(XML_NODE_ROOT);
    }
            
    void Writer4::setTmp(const bool& tmp) const {
        _root.attribute(XML_ATTR_TRANSMEMBRANE).set_value(tmp?"yes":"no");
    }
    
    void Writer4::setCode(const std::string& code) const {
        _root.attribute(XML_ATTR_PDB_CODE).set_value(code.c_str());
    }
            
    void Writer4::setCreateDate(const std::string& date) const {
        _root.child(XML_NODE_CREATED).text() = date.c_str();
    }

    void Writer4::setModifications(const std::vector<Tmdet::ValueObjects::Modification>& mods) {
        pugi::xml_node node;
        for(const auto& m : mods) {
            node = _root.insert_child_after(XML_NODE_MODIFICATION, _root.child(XML_NODE_RAWDATA));
            node.append_attribute(XML_ATTR_DATE);
            node.attribute(XML_ATTR_DATE) = m.date.c_str();
            node.set_value(m.descr.c_str());
        }
    }

    void Writer4::setQvalue(const double& q) const {
        _root.child(XML_NODE_RAWDATA).child(XML_NODE_QVALUE).text() = std::format("{:.2f}",q).c_str();
    }

    void Writer4::setTmtype(const std::string& ptype) const {
        _root.child(XML_NODE_RAWDATA).child(XML_NODE_TMTYPE).text() = ptype.c_str();
    }

    void Writer4::setTMatrix(const pugi::xml_node& node, Tmdet::ValueObjects::TMatrix& tmatrix) const {
        node.child(XML_NODE_TRANSFORMATION)
            .child(XML_NODE_ROTATE)
            .child(XML_NODE_ROWX)
            .attribute(XML_ATTR_X) = std::to_string(tmatrix.rot[0][0]).c_str();
        node.child(XML_NODE_TRANSFORMATION)
            .child(XML_NODE_ROTATE)
            .child(XML_NODE_ROWX)
            .attribute(XML_ATTR_Y) = std::to_string(tmatrix.rot[0][1]).c_str();
        node.child(XML_NODE_TRANSFORMATION)
            .child(XML_NODE_ROTATE)
            .child(XML_NODE_ROWX)
            .attribute(XML_ATTR_Z) = std::to_string(tmatrix.rot[0][2]).c_str();
        node.child(XML_NODE_TRANSFORMATION)
            .child(XML_NODE_ROTATE)
            .child(XML_NODE_ROWY)
            .attribute(XML_ATTR_X) = std::to_string(tmatrix.rot[1][0]).c_str();
        node.child(XML_NODE_TRANSFORMATION)
            .child(XML_NODE_ROTATE)
            .child(XML_NODE_ROWY)
            .attribute(XML_ATTR_Y) = std::to_string(tmatrix.rot[1][1]).c_str();
        node.child(XML_NODE_TRANSFORMATION)
            .child(XML_NODE_ROTATE)
            .child(XML_NODE_ROWY)
            .attribute(XML_ATTR_Z) = std::to_string(tmatrix.rot[1][2]).c_str();
        node.child(XML_NODE_TRANSFORMATION)
            .child(XML_NODE_ROTATE)
            .child(XML_NODE_ROWZ)
            .attribute(XML_ATTR_X) = std::to_string(tmatrix.rot[2][0]).c_str();
        node.child(XML_NODE_TRANSFORMATION)
            .child(XML_NODE_ROTATE)
            .child(XML_NODE_ROWZ)
            .attribute(XML_ATTR_Y) = std::to_string(tmatrix.rot[2][1]).c_str();
        node.child(XML_NODE_TRANSFORMATION)
            .child(XML_NODE_ROTATE)
            .child(XML_NODE_ROWZ)
            .attribute(XML_ATTR_Z) = std::to_string(tmatrix.rot[2][2]).c_str();
        node.child(XML_NODE_TRANSFORMATION)
            .child(XML_NODE_TRANSLATE)
            .attribute(XML_ATTR_X) = std::to_string(tmatrix.trans.x).c_str();
        node.child(XML_NODE_TRANSFORMATION)
            .child(XML_NODE_TRANSLATE)
            .attribute(XML_ATTR_Y) = std::to_string(tmatrix.trans.y).c_str();
        node.child(XML_NODE_TRANSFORMATION)
            .child(XML_NODE_TRANSLATE)
            .attribute(XML_ATTR_Z) = std::to_string(tmatrix.trans.z).c_str();
    }

    void Writer4::setMembranes(std::vector<Tmdet::ValueObjects::Membrane>& membranes) {
        
        for(auto& m : membranes) {
            pugi::xml_document doc;
            doc.load_string(_membrane_xml.c_str());
            pugi::xml_node node = doc.child(XML_NODE_MEMBRANE);
            node.attribute(XML_ATTR_HALF_THICKNESS) = std::to_string(m.halfThickness).c_str();
            node.attribute(XML_ATTR_TYPE) = m.type.name.c_str();
            //setTMatrix(node,m.tmatrix);
            _root.append_copy(node);
        }
    }

    void Writer4::setChains(const std::vector<Tmdet::ValueObjects::Chain>& chains) {
        for(const auto& c: chains) {
            pugi::xml_document doc;
            doc.load_string(_chain_xml.c_str());
            pugi::xml_node node = doc.child(XML_NODE_CHAIN);
            node.attribute(XML_ATTR_AUTH_ASYM_ID) = c.id.c_str();
            node.attribute(XML_ATTR_LABEL_ASYM_ID) = c.labId.c_str();
            node.attribute(XML_ATTR_NUM_TM) = std::to_string(c.numtm).c_str();
            node.attribute(XML_ATTR_TYPE) = c.type.name.c_str();
            node.child(XML_NODE_SEQENCE).text() = c.seq.c_str();
            setRegions(node,c.regions);
            _root.append_copy(node);
        }
    }

    void Writer4::setRegions(pugi::xml_node& pnode, const std::vector<Tmdet::ValueObjects::Region>& regions) const {
        for(const auto& r: regions) {
            pugi::xml_document doc;
            doc.load_string(_region_xml.c_str());
            pugi::xml_node node = doc.child(XML_NODE_REGION);
            
            node.child(XML_NODE_START)
                .attribute(XML_ATTR_SEQ) = std::to_string(r.beg).c_str();
            node.child(XML_NODE_START)
                .attribute(XML_ATTR_AUTH_SEQ_ID) = std::to_string(r.beg).c_str();
            node.child(XML_NODE_START)
                .attribute(XML_ATTR_LABEL_SEQ_ID) = std::to_string(r.beg).c_str();
            node.child(XML_NODE_END)
                .attribute(XML_ATTR_SEQ) = std::to_string(r.end).c_str();
            node.child(XML_NODE_END)
                .attribute(XML_ATTR_AUTH_SEQ_ID) = std::to_string(r.end).c_str();
            node.child(XML_NODE_END)
                .attribute(XML_ATTR_LABEL_SEQ_ID) = std::to_string(r.end).c_str();
            node.attribute(XML_ATTR_TYPE) = std::format("{}","%c",r.type.code).c_str();
            pnode.append_copy(node);
        }
    }
            
    void Writer4::writeXml(Tmdet::ValueObjects::Protein& protein, const std::string& path) {
        create();
        setTmp(protein.tmp);
        setCode(protein.code);
        setCreateDate(protein.date);
        setModifications(protein.modifications);
        setQvalue(protein.qValue);
        setTmtype(protein.type.name);
        setMembranes(protein.membranes);
        setChains(protein.chains);
        write(path);
    }

}
