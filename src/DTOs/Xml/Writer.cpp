#include <string>
#include <vector>
#include <iostream>
#include <format>
#include <pugixml.hpp>
#include <Config.hpp>
#include <System/Logger.hpp>
#include <DTOs/Xml/Writer.hpp>
#include <DTOs/Xml/Constants4.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/TMatrix.hpp>
#include <ValueObjects/Chain.hpp>
#include <Exceptions/SyntaxErrorException.hpp>

namespace Tmdet::DTOs::Xml {

    void Writer::write(const std::string& path) const {
        _doc.save_file(path.c_str(),"  ");
    }

    void Writer::create() {
        _doc.load_string(_pdbtm_xml.c_str());
        _root = _doc.child(XML_NODE_ROOT);
    }
            
    void Writer::setTmp(const bool& tmp) const {
        _root.attribute(XML_ATTR_TRANSMEMBRANE).set_value(tmp?"yes":"no");
    }
    
    void Writer::setCode(const std::string& code) const {
        _root.attribute(XML_ATTR_PDB_CODE).set_value(code.c_str());
    }
            
    void Writer::setCreateDate(const std::string& date) const {
        _root.child(XML_NODE_RAWDATA).child(XML_NODE_CREATED).text() = date.c_str();
    }

    void Writer::setVersion(const std::string& version) const {
        _root.child(XML_NODE_RAWDATA).child(XML_NODE_TMDET_VERSION).text() = version.c_str();
    }

    void Writer::setQvalue(const double& q) const {
        _root.child(XML_NODE_RAWDATA).child(XML_NODE_QVALUE).text() = std::format("{:.2f}",q).c_str();
    }

    void Writer::setTmtype(const std::string& ptype) const {
        _root.child(XML_NODE_RAWDATA).child(XML_NODE_TMTYPE).text() = ptype.c_str();
    }

    void Writer::setTMatrix(Tmdet::ValueObjects::TMatrix& tmatrix) {
        auto pnode = _root.insert_child_after(XML_NODE_TRANSFORMATION, _root.child(XML_NODE_RAWDATA));
        pugi::xml_node node = pnode.append_child(XML_NODE_TRANSLATE);
        node.append_attribute(XML_ATTR_X) = std::format("{:.6f}",tmatrix.trans.x).c_str();
        node.append_attribute(XML_ATTR_Y) = std::format("{:.6f}",tmatrix.trans.y).c_str();
        node.append_attribute(XML_ATTR_Z) = std::format("{:.6f}",tmatrix.trans.z).c_str();
        node = pnode.append_child(XML_NODE_ROTATE);
        pugi::xml_node rowNode = node.append_child(XML_NODE_ROWX);
        rowNode.append_attribute(XML_ATTR_X) = std::format("{:.6f}",tmatrix.rot[0][0]).c_str();
        rowNode.append_attribute(XML_ATTR_Y) = std::format("{:.6f}",tmatrix.rot[0][1]).c_str();
        rowNode.append_attribute(XML_ATTR_Z) = std::format("{:.6f}",tmatrix.rot[0][2]).c_str();
        rowNode = node.append_child(XML_NODE_ROWY);
        rowNode.append_attribute(XML_ATTR_X) = std::format("{:.6f}",tmatrix.rot[1][0]).c_str();
        rowNode.append_attribute(XML_ATTR_Y) = std::format("{:.6f}",tmatrix.rot[1][1]).c_str();
        rowNode.append_attribute(XML_ATTR_Z) = std::format("{:.6f}",tmatrix.rot[1][2]).c_str();
        rowNode = node.append_child(XML_NODE_ROWZ);
        rowNode.append_attribute(XML_ATTR_X) = std::format("{:.6f}",tmatrix.rot[2][0]).c_str();
        rowNode.append_attribute(XML_ATTR_Y) = std::format("{:.6f}",tmatrix.rot[2][1]).c_str();
        rowNode.append_attribute(XML_ATTR_Z) = std::format("{:.6f}",tmatrix.rot[2][2]).c_str();
    }

    void Writer::setMembranes(std::vector<Tmdet::ValueObjects::Membrane>& membranes) {
        auto pnode = _root.insert_child_after(XML_NODE_MEMBRANES, _root.child(XML_NODE_TRANSFORMATION));
        for(auto& membrane : membranes) {
            pugi::xml_node node = pnode.append_child(XML_NODE_MEMBRANE);
            node.append_attribute(XML_ATTR_HALF_THICKNESS) = std::format("{:.1f}",membrane.halfThickness).c_str();
            node.append_attribute(XML_ATTR_SIZE) = std::format("{:.1f}",membrane.membraneRadius).c_str();
            node.append_attribute(XML_ATTR_TYPE) = membrane.type.name.c_str();
            node.append_attribute(XML_ATTR_Z) = std::format("{:.1f}",membrane.origo).c_str();
        }
    }

    void Writer::setChains(const std::vector<Tmdet::ValueObjects::Chain>& chains) {
        auto pnode = _root.insert_child_after(XML_NODE_CHAINS, _root.child(XML_NODE_MEMBRANES));
        for(const auto& chain: chains) {
            pugi::xml_node node = pnode.append_child(XML_NODE_CHAIN);
            node.append_attribute(XML_ATTR_AUTH_ASYM_ID) = chain.id.c_str();
            node.append_attribute(XML_ATTR_LABEL_ASYM_ID) = chain.labId.c_str();
            node.append_attribute(XML_ATTR_NUM_TM) = std::to_string(chain.numtm).c_str();
            node.append_attribute(XML_ATTR_TYPE) = chain.type.name.c_str();
            pugi::xml_node seqNode = node.append_child(XML_NODE_SEQENCE);
            seqNode.text() = chain.seq.c_str(); //TODO format sequence
            if (chain.selected && !chain.regions.empty()) {
                setRegions(node,chain.regions);
            }
        }
    }

    void Writer::setRegions(pugi::xml_node& pnode, const std::vector<Tmdet::ValueObjects::Region>& regions) const {
        for(const auto& r: regions) {
            pugi::xml_node node = pnode.append_child(XML_NODE_REGION);
            pugi::xml_node startNode = node.append_child(XML_NODE_START);
            pugi::xml_node endNode = node.append_child(XML_NODE_END);
            startNode.append_attribute(XML_ATTR_SEQ) = std::to_string(r.beg).c_str();
            startNode.append_attribute(XML_ATTR_AUTH_SEQ_ID) = std::to_string(r.beg_auth_seq_id).c_str();
            if (r.beg_auth_seq_icode != ' ') {
                startNode.append_attribute(XML_ATTR_AUTH_SEQ_ICODE) = std::to_string(r.beg_auth_seq_icode).c_str();
            }
            startNode.append_attribute(XML_ATTR_LABEL_SEQ_ID) = std::to_string(r.beg_label_seq_id).c_str();
            endNode.append_attribute(XML_ATTR_SEQ) = std::to_string(r.end).c_str();
            endNode.append_attribute(XML_ATTR_AUTH_SEQ_ID) = std::to_string(r.end_auth_seq_id).c_str();
            if (r.end_auth_seq_icode !=  ' ')  {
                endNode.append_attribute(XML_ATTR_AUTH_SEQ_ID) = std::to_string(r.end_auth_seq_icode).c_str();
            }
            endNode.append_attribute(XML_ATTR_LABEL_SEQ_ID) = std::to_string(r.end_label_seq_id).c_str();

            node.append_attribute(XML_ATTR_TYPE) = std::format("{}",r.type.code).c_str();
        }
    }
            
    void Writer::writeXml(Tmdet::ValueObjects::Protein& protein, const std::string& path) {
        DEBUG_LOG("Processing: Writer::writeXml({})",path);
        create();
        setTmp(protein.tmp);
        setCode(protein.code);
        setCreateDate(protein.date);
        setVersion(protein.version);
        setQvalue(protein.qValue);
        setTmtype(protein.type.name);
        setTMatrix(protein.tmatrix);
        setMembranes(protein.membranes);
        setChains(protein.chains);
        write(path);
        DEBUG_LOG(" Processed: Writer::writeXml({})",path);
    }

}
