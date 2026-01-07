// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <string>
#include <vector>
#include <iostream>
#include <format>
#include <pugixml.hpp>
#include <Config.hpp>
#include <DTOs/XmlRW/Constants4.hpp>
#include <DTOs/XmlRW/Writer.hpp>
#include <Exceptions/SyntaxErrorException.hpp>
#include <Helpers/String.hpp>
#include <System/Arguments.hpp>
#include <System/Logger.hpp>
#include <VOs/Membrane.hpp>
#include <VOs/TMatrix.hpp>
#include <VOs/Chain.hpp>
#include <VOs/Xml.hpp>

namespace Tmdet::DTOs::XmlRW {

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
        DEBUG_LOG("setCode: >>{}<<",code);
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

    void Writer::setArguments(const Tmdet::System::Arguments& args) const {
        _root.child(XML_NODE_RAWDATA).child(XML_NODE_ARGUMENTS).text() = args.getCommandLine().c_str();
    }

    void Writer::setTMatrix(Tmdet::VOs::TMatrix& tmatrix) {
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

    void Writer::setMembranes(std::vector<Tmdet::VOs::Membrane>& membranes) {
        auto pnode = _root.insert_child_after(XML_NODE_MEMBRANES, _root.child(XML_NODE_TRANSFORMATION));
        for(auto& membrane : membranes) {
            pugi::xml_node node = pnode.append_child(XML_NODE_MEMBRANE);
            node.append_attribute(XML_ATTR_HALF_THICKNESS) = std::format("{:.1f}",membrane.halfThickness).c_str();
            node.append_attribute(XML_ATTR_SIZE) = std::format("{:.1f}",membrane.membraneRadius).c_str();
            node.append_attribute(XML_ATTR_TYPE) = membrane.type.name.c_str();
            node.append_attribute(XML_ATTR_Z) = std::format("{:.1f}",membrane.origo).c_str();
            if (!membrane.type.isPlane()) {
                node.append_attribute(XML_ATTR_SPHERE_RADIUS) = std::format("{:.1f}",membrane.sphereRadius).c_str();
            }
        }
    }

    void Writer::setChains(const std::vector<Tmdet::VOs::XmlChain>& chains) {
        auto pnode = _root.insert_child_after(XML_NODE_CHAINS, _root.child(XML_NODE_MEMBRANES));
        for(const auto& chain: chains) {
            pugi::xml_node node = pnode.append_child(XML_NODE_CHAIN);
            node.append_attribute(XML_ATTR_AUTH_ID) = chain.id.c_str();
            node.append_attribute(XML_ATTR_LABEL_ID) = chain.labId.c_str();
            node.append_attribute(XML_ATTR_NUM_TM) = std::to_string(chain.numtm).c_str();
            node.append_attribute(XML_ATTR_TYPE) = chain.type.name.c_str();
            pugi::xml_node seqNode = node.append_child(XML_NODE_SEQENCE);
            seqNode.text() = ((std::string)"\n"+Tmdet::Helpers::String::formatSequence(chain.seq,50,10,"        ")+(std::string)"\n      ").c_str();
            if (chain.selected && !chain.regions.empty()) {
                setRegions(node, chain.regions);
            }
        }
    }

    void Writer::setRegions(pugi::xml_node& pnode, const std::vector<Tmdet::VOs::Region>& regions) const {
        pugi::xml_node regions_node = pnode.append_child(XML_NODE_REGIONS);
        for(const auto& r: regions) {
            pugi::xml_node node = regions_node.append_child(XML_NODE_REGION);
            node.append_attribute(XML_ATTR_START_AUTH_ID) = std::to_string(r.beg.authId).c_str();
            if (r.beg.authIcode != ' ') {
                node.append_attribute(XML_ATTR_START_AUTH_ICODE) = std::to_string(r.beg.authIcode).c_str();
            }
            node.append_attribute(XML_ATTR_START_LABEL_ID) = std::to_string(r.beg.labelId).c_str();
            node.append_attribute(XML_ATTR_END_AUTH_ID) = std::to_string(r.end.authId).c_str();
            if (r.end.authIcode !=  ' ')  {
                node.append_attribute(XML_ATTR_END_AUTH_ICODE) = std::to_string(r.end.authIcode).c_str();
            }
            node.append_attribute(XML_ATTR_END_LABEL_ID) = std::to_string(r.end.labelId).c_str();
            node.append_attribute(XML_ATTR_TYPE) = std::format("{}",r.type.code).c_str();
        }
    }
            
    void Writer::writeXml(Tmdet::VOs::Xml& xmlData, const std::string& path, const Tmdet::System::Arguments& args) {
        DEBUG_LOG("Processing: Writer::writeXml({})",path);
        create();
        setTmp(xmlData.tmp);
        setCode(xmlData.code);
        setCreateDate(xmlData.date);
        setVersion(xmlData.version);
        setQvalue(xmlData.qValue);
        setTmtype(xmlData.type.name);
        setArguments(args);
        if (xmlData.tmp) {
            setTMatrix(xmlData.tmatrix);
            setMembranes(xmlData.membranes);
            setChains(xmlData.chains);
        }
        write(path);
        DEBUG_LOG(" Processed: Writer::writeXml({} {})",path,args.getCommandLine());
    }

}
