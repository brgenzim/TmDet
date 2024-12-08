#include <string>
#include <vector>
#include <iostream>
#include <format>
#include <pugixml.hpp>
#include <Config.hpp>
#include <System/Logger.hpp>
#include <DTOs/XmlRW/Constants3.hpp>
#include <DTOs/XmlRW/Reader3.hpp>
#include <VOs/Membrane.hpp>
#include <VOs/BioMatrix.hpp>
#include <VOs/Modification.hpp>
#include <VOs/TMatrix.hpp>
#include <VOs/Xml.hpp>

namespace Tmdet::DTOs::XmlRW {

    void Reader3::setRoot(const pugi::xml_document& doc) {
        _root = doc.child(XML3_NODE_ROOT);
    }

    void Reader3::readXml(Tmdet::VOs::Xml& xmlData) {
        xmlData.tmp = getTmp();
        DEBUG_LOG("tmp: {}",(xmlData.tmp?"yes":"no"));
        xmlData.code = getCode();
        xmlData.date = getCreateDate();
        xmlData.modifications = getModifications();
        xmlData.qValue = getQvalue();
        DEBUG_LOG("qValue: {}",xmlData.qValue);
        auto type = getTmtype();
        DEBUG_LOG("type: {}",type);
        xmlData.type = Tmdet::Types::Proteins.at(type);
        xmlData.spres = getSpres();
        xmlData.pdbkwres = getPdbkwres();
        xmlData.bioMatrix = getBioMatrix();
        xmlData.membranes = getMembranes();
        xmlData.tmatrix = getTMatrix(_root.child(XML3_NODE_MEMBRANE).child(XML3_NODE_TMATRIX));
        xmlData.chains = getChains();
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

    std::vector<Tmdet::VOs::Modification> Reader3::getModifications() const {
        std::vector<Tmdet::VOs::Modification> mods;
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

    Tmdet::VOs::TMatrix Reader3::getTMatrix(const pugi::xml_node& node) const {
        Tmdet::VOs::TMatrix tmatrix;
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

    Tmdet::VOs::BioMatrix Reader3::getBioMatrix() const {
        Tmdet::VOs::BioMatrix bioMatrix;
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

    std::vector<Tmdet::VOs::Membrane> Reader3::getMembranes() const {
        std::vector<Tmdet::VOs::Membrane> membranes;
        for (pugi::xml_node m_node = _root.child(XML3_NODE_MEMBRANE); m_node; m_node = m_node.next_sibling(XML3_NODE_MEMBRANE)) {
            membranes.emplace_back(
                0.0, //todo get origo
                m_node.child(XML3_NODE_NORMAL).attribute(XML3_ATTR_Z).as_double(),
                0.0,
                10.0,
                Tmdet::Types::Membranes.at("Plain")
            );
        }
        return membranes;
    }

    std::vector<Tmdet::VOs::XmlChain> Reader3::getChains() {
        std::vector<Tmdet::VOs::XmlChain> xmlChains;
        for (pugi::xml_node chainNode: _root.children(XML3_NODE_CHAIN)) {
            auto type = chainNode.attribute(XML3_ATTR_TYPE).value();
            xmlChains.emplace_back(
                chainNode.attribute(XML3_ATTR_CHAINID).value(), //id
                chainNode.attribute(XML3_ATTR_CHAINID).value(), //labId
                true, //selected
                chainNode.attribute(XML3_ATTR_NUM_TM).as_int(), //numtm
                chainNode.child(XML3_NODE_SEQ).text().get(), //seq
                getRegions(chainNode), //regions
                Tmdet::Types::Chains.at(type) //type
            );
        }
        return xmlChains;
    }

    std::vector<Tmdet::VOs::Region> Reader3::getRegions(const pugi::xml_node& cnode) const {
        std::vector<Tmdet::VOs::Region> regions;
        for (pugi::xml_node r_node = cnode.child(XML3_NODE_REGION); r_node; r_node = r_node.next_sibling(XML3_NODE_REGION)) {
            char type = r_node.attribute(XML3_ATTR_type).value()[0];
            if (type=='I') {
                type = 'N';
            }
            Tmdet::VOs::Region region = {
                {r_node.attribute(XML3_ATTR_SEQ_BEG).as_int(),' ',r_node.attribute(XML3_ATTR_PDB_BEG).as_int()},
                {r_node.attribute(XML3_ATTR_SEQ_END).as_int(),' ',r_node.attribute(XML3_ATTR_PDB_END).as_int()},
                Tmdet::Types::Regions.at(type)
            };
            regions.push_back(region);
        }
        return regions;
    }

}
