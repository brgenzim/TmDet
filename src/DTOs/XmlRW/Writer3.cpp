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
#include <DTOs/XmlRW/Constants3.hpp>
#include <DTOs/XmlRW/Writer3.hpp>
#include <Exceptions/SyntaxErrorException.hpp>
#include <Helpers/String.hpp>
#include <System/Arguments.hpp>
#include <System/Logger.hpp>
#include <VOs/TMatrix.hpp>
#include <VOs/BioMatrix.hpp>
#include <VOs/Chain.hpp>
#include <VOs/Xml.hpp>

namespace Tmdet::DTOs::XmlRW {

    void Writer3::write(const std::string& path) const {
        _doc.save_file(path.c_str(),"  ");
    }

    void Writer3::create() {
        _doc.load_string(_pdbtm_xml.c_str());
        _root = _doc.child(XML3_NODE_ROOT);
    }
            
    void Writer3::setTmp(const bool& tmp) const {
        _root.attribute(XML3_ATTR_TMP).set_value(tmp?"yes":"no");
    }
    
    void Writer3::setCode(const std::string& code) const {
        _root.attribute(XML3_ATTR_ID).set_value(code.c_str());
    }
            
    void Writer3::setCreateDate(const std::string& date) const {
        _root.child(XML3_NODE_CREATE_DATE).text() = date.c_str();
    }

    void Writer3::setModifications(const std::vector<Tmdet::VOs::Modification>& modifications) {
        for (const auto& mod: modifications) {
            auto node = _root.insert_child_after(XML3_NODE_MODIFICATION, _root.child(XML3_NODE_CREATE_DATE));
            node.append_child(XML3_NODE_DATE).set_value(mod.date.c_str());
            node.append_child(XML3_NODE_DESCR).set_value(mod.descr.c_str()); 
        }
    }

    void Writer3::setQvalue(const double& q) const {
        _root.child(XML3_NODE_RAWRES).child(XML3_NODE_TMRES).text() = std::format("{:.2f}",q).c_str();
    }

    void Writer3::setTmtype(const std::string& ptype) const {
        _root.child(XML3_NODE_RAWRES).child(XML3_NODE_TMTYPE).text() = ptype.c_str();
    }

    void Writer3::setSpres(const std::string& spres) const {
        _root.child(XML3_NODE_RAWRES).child(XML3_NODE_SPRES).text() = spres.c_str();
    }

    void Writer3::setPdbkw(const std::string& pdbkw) const {
        _root.child(XML3_NODE_RAWRES).child(XML3_NODE_PDBKWRES).text() = pdbkw.c_str();
    }

    void Writer3::setArguments(const Tmdet::System::Arguments& args) const {
        
    }

    void Writer3::setTMatrix(pugi::xml_node& pnode, Tmdet::VOs::TMatrix& tmatrix) {
        
        pugi::xml_node node = pnode.append_child(XML3_NODE_ROWX);
        node.append_attribute(XML3_ATTR_X) = std::format("{:.6f}",tmatrix.rot[0][0]).c_str();
        node.append_attribute(XML3_ATTR_Y) = std::format("{:.6f}",tmatrix.rot[0][1]).c_str();
        node.append_attribute(XML3_ATTR_Z) = std::format("{:.6f}",tmatrix.rot[0][2]).c_str();
        node.append_attribute(XML3_ATTR_T) = std::format("{:.6f}",tmatrix.trans.x).c_str();
        node = pnode.append_child(XML3_NODE_ROWY);
        node.append_attribute(XML3_ATTR_X) = std::format("{:.6f}",tmatrix.rot[1][0]).c_str();
        node.append_attribute(XML3_ATTR_Y) = std::format("{:.6f}",tmatrix.rot[1][1]).c_str();
        node.append_attribute(XML3_ATTR_Z) = std::format("{:.6f}",tmatrix.rot[1][2]).c_str();
        node.append_attribute(XML3_ATTR_T) = std::format("{:.6f}",tmatrix.trans.y).c_str();
        node = pnode.append_child(XML3_NODE_ROWZ);
        node.append_attribute(XML3_ATTR_X) = std::format("{:.6f}",tmatrix.rot[2][0]).c_str();
        node.append_attribute(XML3_ATTR_Y) = std::format("{:.6f}",tmatrix.rot[2][1]).c_str();
        node.append_attribute(XML3_ATTR_Z) = std::format("{:.6f}",tmatrix.rot[2][2]).c_str();
        node.append_attribute(XML3_ATTR_T) = std::format("{:.6f}",tmatrix.trans.x).c_str();
    }

    void Writer3::setBioMatrix(const Tmdet::VOs::BioMatrix& bioMatrix) {

    }

    void Writer3::setMembrane(std::vector<Tmdet::VOs::Membrane>& membranes, Tmdet::VOs::TMatrix& tmatrix) {
        auto pnode = _root.insert_child_after(XML3_NODE_MEMBRANE,(_root.child(XML3_NODE_BIOMATRIX)?_root.child(XML3_NODE_BIOMATRIX):_root.child(XML3_NODE_RAWRES)));
        auto node = pnode.append_child(XML3_NODE_NORMAL);
        node.append_attribute(XML3_ATTR_X) = "0.0000";
        node.append_attribute(XML3_ATTR_Y) = "0.0000";
        node.append_attribute(XML3_ATTR_Z) = std::format("{:.6f}",membranes[0].halfThickness).c_str();
        node = pnode.append_child(XML3_NODE_TMATRIX);
        setTMatrix(node,tmatrix);
    }

    void Writer3::setChains(const std::vector<Tmdet::VOs::XmlChain>& chains) {
        for(const auto& chain: chains) {
            pugi::xml_node node = _root.insert_child_after(XML3_NODE_CHAIN, (_root.child(XML3_NODE_CHAIN)?_root.child(XML3_NODE_CHAIN):_root.child(XML3_NODE_MEMBRANE)));
            node.append_attribute(XML3_ATTR_CHAINID) = chain.id.c_str();
            node.append_attribute(XML3_ATTR_NUM_TM) = std::to_string(chain.numtm).c_str();
            node.append_attribute(XML3_ATTR_TYPE) = chain.type.name.c_str();
            pugi::xml_node seqNode = node.append_child(XML3_NODE_SEQ);
            seqNode.text() = ((std::string)"\n"+Tmdet::Helpers::String::formatSequence(chain.seq,50,10,"        ")+(std::string)"\n      ").c_str();
            if (chain.selected && !chain.regions.empty()) {
                setRegions(node, chain.regions);
            }
        }
    }

    void Writer3::setRegions(pugi::xml_node& pnode, const std::vector<Tmdet::VOs::Region>& regions) const {
        for(const auto& r: regions) {
            pugi::xml_node node = pnode.append_child(XML3_NODE_REGION);
            node.append_attribute(XML3_ATTR_SEQ_BEG) = std::to_string(r.beg.labelId).c_str();
            node.append_attribute(XML3_ATTR_PDB_BEG) = std::to_string(r.beg.authId).c_str();
            if (r.beg.authIcode != ' ') {
                node.append_attribute(XML3_ATTR_PDB_BEGI) = std::to_string(r.beg.authIcode).c_str();
            }
            node.append_attribute(XML3_ATTR_SEQ_END) = std::to_string(r.end.labelId).c_str();
            node.append_attribute(XML3_ATTR_PDB_END) = std::to_string(r.end.authId).c_str();
            if (r.end.authIcode != ' ') {
                node.append_attribute(XML3_ATTR_PDB_ENDI) = std::to_string(r.end.authIcode).c_str();
            }
            node.append_attribute(XML3_ATTR_type) = std::format("{}",r.type.code).c_str();
        }
    }
            
    void Writer3::writeXml(Tmdet::VOs::Xml& xmlData, const std::string& path, const Tmdet::System::Arguments& args) {
        DEBUG_LOG("Processing: Writer3::writeXml({})",path);
        create();
        setTmp(xmlData.tmp);
        setCode(xmlData.code);
        setCreateDate(xmlData.date);
        setModifications(xmlData.modifications);
        setQvalue(xmlData.qValue);
        setTmtype(xmlData.type.name);
        setSpres(xmlData.spres);
        setPdbkw(xmlData.pdbkwres);
        setArguments(args);
        if (xmlData.tmp) {
            setBioMatrix(xmlData.bioMatrix);
            setMembrane(xmlData.membranes, xmlData.tmatrix);
            setChains(xmlData.chains);
        }
        write(path);
        DEBUG_LOG(" Processed: Writer3::writeXml({})",path);
    }

}
