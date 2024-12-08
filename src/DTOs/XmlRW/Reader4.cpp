#include <string>
#include <vector>
#include <pugixml.hpp>
#include <Config.hpp>
#include <DTOs/Region.hpp>
#include <DTOs/XmlRW/Constants4.hpp>
#include <DTOs/XmlRW/Reader4.hpp>
#include <System/Logger.hpp>
#include <VOs/Membrane.hpp>
#include <VOs/BioMatrix.hpp>
#include <VOs/Modification.hpp>
#include <VOs/TMatrix.hpp>
#include <VOs/Xml.hpp>

namespace Tmdet::DTOs::XmlRW {

    void Reader4::setRoot(const pugi::xml_document& doc) {
        _root = doc.child(XML_NODE_ROOT);
    }

    void Reader4::readXml(Tmdet::VOs::Xml& xmlData) {
        xmlData.tmp = getTmp();
        xmlData.code = getCode();
        xmlData.date = getCreateDate();
        xmlData.version = getVersion();
        xmlData.qValue = getQvalue();
        xmlData.type = Tmdet::Types::Proteins.at(getTmtype());
        xmlData.membranes = getMembranes();
        xmlData.chains = getChains();
    }

    bool Reader4::getTmp() const {
        return _root.attribute(XML_ATTR_TRANSMEMBRANE).as_bool();
    }

    std::string Reader4::getCode() const {
        return _root.attribute(XML_ATTR_PDB_CODE).value();
    }

    std::string Reader4::getCreateDate() const {
        return _root.child(XML_NODE_RAWDATA).child(XML_NODE_CREATED).text().get();
    }

    std::string Reader4::getVersion() const {
        return _root.child(XML_NODE_RAWDATA).child(XML_NODE_TMDET_VERSION).text().get();
    }

    double Reader4::getQvalue() const {
        return _root.child(XML_NODE_RAWDATA).child(XML_NODE_QVALUE).text().as_double();
    }

    std::string Reader4::getTmtype() const {
        return _root.child(XML_NODE_RAWDATA).child(XML_NODE_TMTYPE).text().get();
    }

    Tmdet::VOs::TMatrix Reader4::getTMatrix(const pugi::xml_node& pnode) const {
        Tmdet::VOs::TMatrix tmatrix;
        pugi::xml_node node = pnode.child(XML_NODE_TRANSLATE);
        tmatrix.trans.x = node.attribute(XML_ATTR_X).as_double();
        tmatrix.trans.y = node.attribute(XML_ATTR_Y).as_double();
        tmatrix.trans.z = node.attribute(XML_ATTR_Z).as_double();
        node = pnode.child(XML_NODE_ROTATE);
        pugi::xml_node rowNode = node.child(XML_NODE_ROWX);
        tmatrix.rot[0][0] = rowNode.attribute(XML_ATTR_X).as_double();
        tmatrix.rot[0][1] = rowNode.attribute(XML_ATTR_Y).as_double();
        tmatrix.rot[0][2] = rowNode.attribute(XML_ATTR_Z).as_double();
        rowNode = node.child(XML_NODE_ROWY);
        tmatrix.rot[1][0] = rowNode.attribute(XML_ATTR_X).as_double();
        tmatrix.rot[1][1] = rowNode.attribute(XML_ATTR_Y).as_double();
        tmatrix.rot[1][2] = rowNode.attribute(XML_ATTR_Z).as_double();
        rowNode = node.child(XML_NODE_ROWZ);
        tmatrix.rot[2][0] = rowNode.attribute(XML_ATTR_X).as_double();
        tmatrix.rot[2][1] = rowNode.attribute(XML_ATTR_Y).as_double();
        tmatrix.rot[2][2] = rowNode.attribute(XML_ATTR_Z).as_double();
        return tmatrix;
    }

    std::vector<Tmdet::VOs::Membrane> Reader4::getMembranes() const {
        auto pnode = _root.child(XML_NODE_MEMBRANES);
        std::vector<Tmdet::VOs::Membrane> membranes;
        for (pugi::xml_node membraneNode: pnode.children(XML_NODE_MEMBRANE)) {
            membranes.emplace_back(
                membraneNode.attribute(XML_ATTR_Z).as_double(),
                membraneNode.attribute(XML_ATTR_HALF_THICKNESS).as_double(),
                membraneNode.attribute(XML_ATTR_SPHERE_RADIUS)?
                    membraneNode.attribute(XML_ATTR_SPHERE_RADIUS).as_double():0.0),
                membraneNode.attribute(XML_ATTR_SIZE).as_double(),
                Tmdet::Types::Membranes.at(membraneNode.attribute(XML_ATTR_TYPE).as_string()
            );
        }
        return membranes;
    }

    std::vector<Tmdet::VOs::XmlChain> Reader4::getChains() {
        std::vector<Tmdet::VOs::XmlChain> xmlChains;
        auto pnode = _root.child(XML_NODE_CHAINS);
        for (pugi::xml_node chainNode: pnode.children(XML_NODE_CHAIN)) {
            auto type = chainNode.attribute(XML_ATTR_TYPE).as_string();
            xmlChains.emplace_back(
                chainNode.attribute(XML_ATTR_AUTH_ID).as_string(),
                chainNode.attribute(XML_ATTR_LABEL_ID).as_string(),
                (type != Tmdet::Types::ChainType::NOT_SELECTED.name),
                chainNode.attribute(XML_ATTR_NUM_TM).as_int(),
                chainNode.child(XML_NODE_SEQENCE).text().get(),
                getRegions(chainNode),
                Tmdet::Types::Chains.at(type)
            );
        }
        return xmlChains;
    }

    std::vector<Tmdet::VOs::Region> Reader4::getRegions(const pugi::xml_node& pnode) const {
        std::vector<Tmdet::VOs::Region> regions;
        
        for(pugi::xml_node node: pnode.child(XML_NODE_REGIONS).children(XML_NODE_REGION)) {
            Tmdet::VOs::Region region = {
                {
                    node.attribute(XML_ATTR_START_AUTH_ID).as_int(),
                    (node.attribute(XML_ATTR_START_AUTH_ICODE)?node.attribute(XML_ATTR_START_AUTH_ICODE).as_string()[0]:' '),
                    node.attribute(XML_ATTR_START_LABEL_ID).as_int()
                },
                {
                    node.attribute(XML_ATTR_END_AUTH_ID).as_int(),
                    (node.attribute(XML_ATTR_END_AUTH_ICODE)?node.attribute(XML_ATTR_END_AUTH_ICODE).as_string()[0]:' '),
                    node.attribute(XML_ATTR_END_LABEL_ID).as_int()
                },
                Tmdet::Types::Regions.at(node.attribute(XML_ATTR_TYPE).as_string()[0])
            };
            regions.push_back(region);
        }
        return regions;
    }
}