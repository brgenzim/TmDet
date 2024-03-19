#include <string>
#include <vector>
#include <iostream>
#include <pugixml.hpp>
#include <Utils/Xml.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/TMatrix.hpp>
#include <ValueObjects/Chain.hpp>

using namespace std;
using namespace gemmi;

namespace Tmdet::Utils {
    
    void Xml::read(string& path) {
        pugi::xml_parse_result result = _doc.load_file(path.c_str());
        if (!result) {
            cerr << "Error description: " << result.description() << endl;
            cerr << "Error offset: " << result.offset << endl;
            exit(EXIT_FAILURE);
        }
        _root = _doc.child(XML_NODE_ROOT);
    }

    void Xml::write(string& path) {
        _doc.save_file(path.c_str(),"  ");
    }

    void Xml::create() {
        _doc.load_string(_pdbtm_xml.c_str());
        _root = _doc.child(XML_NODE_ROOT);
    }

    bool Xml::getTmp() {
        return _root.attribute(XML_ATTR_TMP).as_bool();
    }

    void Xml::setTmp(bool& tmp) {
        _root.attribute(XML_ATTR_TMP).set_value(tmp?"yes":"no");
    }

    string Xml::getCode() {
        return _root.attribute(XML_ATTR_ID).value();
    }

    void Xml::setCode(string& code) {
        _root.attribute(XML_ATTR_ID).set_value(code.c_str());
    }

    string Xml::getCreateDate() {
        return _root.child(XML_NODE_CREATE_DATE).text().get();
    }

    void Xml::setCreateDate(string& date) {
        _root.child(XML_NODE_CREATE_DATE).text() = date.c_str();
    }

    vector<Tmdet::ValueObjects::Modification> Xml::getModifications() {
        vector<Tmdet::ValueObjects::Modification> mods;
        for (pugi::xml_node mod = _root.child(XML_NODE_MODIFICATION); mod; mod = mod.next_sibling(XML_NODE_MODIFICATION)) {
            Tmdet::ValueObjects::Modification m = {mod.child(XML_NODE_DATE).text().get(),mod.child(XML_NODE_DESCR).text().get()};
            mods.emplace_back(m);
        }
        return mods;
    }

    void Xml::setModifications(vector<Tmdet::ValueObjects::Modification>& mods) {
        pugi::xml_node node;
        for(auto m : mods) {
            node = _root.insert_child_after(XML_NODE_MODIFICATION, _root.child(XML_NODE_CREATE_DATE));
            pugi::xml_node date_node = node.append_child(XML_NODE_DATE);
            date_node.append_child(pugi::node_pcdata).set_value(m.date.c_str());
            pugi::xml_node descr_node = node.append_child(XML_NODE_DESCR);
            descr_node.append_child(pugi::node_pcdata).set_value(m.descr.c_str());
        }
    }

    double Xml::getQvalue() {
        return _root.child(XML_NODE_RAWRES).child(XML_NODE_TMRES).text().as_double();
    }

    void Xml::setQvalue(double& q) {
        char s[100];
        sprintf(s,"%.2f",q);
        _root.child(XML_NODE_RAWRES).child(XML_NODE_TMRES).text() = s;
    }

    string Xml::getTmtype() {
        return _root.child(XML_NODE_RAWRES).child(XML_NODE_TMTYPE).text().get();
    }

    void Xml::setTmtype(string& ptype) {
        _root.child(XML_NODE_RAWRES).child(XML_NODE_TMTYPE).text() = ptype.c_str();
    }

    string Xml::getSpres() {
        return _root.child(XML_NODE_RAWRES).child(XML_NODE_SPRES).text().get();
    }

    void Xml::setSpres(string& spres) {
        _root.child(XML_NODE_RAWRES).child(XML_NODE_SPRES).text() = spres.c_str();
    }

    string Xml::getPdbkwres() {
        return _root.child(XML_NODE_RAWRES).child(XML_NODE_PDBKWRES).text().get();
    }

    void Xml::setPdbkwres(string& pdbkwres) {
        _root.child(XML_NODE_RAWRES).child(XML_NODE_PDBKWRES).text() = pdbkwres.c_str();
    }

    Tmdet::ValueObjects::TMatrix Xml::getTMatrix(pugi::xml_node& node) {
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

    void Xml::setTMatrix(pugi::xml_node& node, Tmdet::ValueObjects::TMatrix& tmatrix) {
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWX).attribute(XML_ATTR_X) = to_string(tmatrix.rot[0][0]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWX).attribute(XML_ATTR_Y) = to_string(tmatrix.rot[0][1]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWX).attribute(XML_ATTR_Z) = to_string(tmatrix.rot[0][2]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWY).attribute(XML_ATTR_X) = to_string(tmatrix.rot[1][0]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWY).attribute(XML_ATTR_Y) = to_string(tmatrix.rot[1][1]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWY).attribute(XML_ATTR_Z) = to_string(tmatrix.rot[1][2]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWZ).attribute(XML_ATTR_X) = to_string(tmatrix.rot[2][0]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWZ).attribute(XML_ATTR_Y) = to_string(tmatrix.rot[2][1]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWZ).attribute(XML_ATTR_Z) = to_string(tmatrix.rot[2][2]).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWX).attribute(XML_ATTR_T) = to_string(tmatrix.trans.x).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWY).attribute(XML_ATTR_T) = to_string(tmatrix.trans.y).c_str();
        node.child(XML_NODE_TMATRIX).child(XML_NODE_ROWZ).attribute(XML_ATTR_T) = to_string(tmatrix.trans.z).c_str();
    }

     Tmdet::ValueObjects::BioMatrix Xml::getBioMatrix() {
        Tmdet::ValueObjects::BioMatrix bioMatrix;
        pugi::xml_node node = _root.child(XML_NODE_BIOMATRIX);
        for (pugi::xml_node matrix = node.child(XML_NODE_MATRIX); matrix; matrix = matrix.next_sibling(XML_NODE_MATRIX)) {
            pugi::xml_node tnode = matrix.child(XML_NODE_TMATRIX);
            Tmdet::ValueObjects::Matrix m = {
                matrix.attribute(XML_ATTR_ID).as_int(),
                matrix.child(XML_NODE_APPLY_TO_CHAIN).attribute(XML_ATTR_CHAINID).as_string(),
                matrix.child(XML_NODE_APPLY_TO_CHAIN).attribute(XML_ATTR_NEW_CHAINID).as_string(),
                getTMatrix(tnode)
            };
            bioMatrix.matrices.emplace_back(m);
        }
        return bioMatrix;
     }

    void Xml::setBioMatrix(Tmdet::ValueObjects::BioMatrix& bioMatrix) {
        pugi::xml_node biom_node = _root.insert_child_after(XML_NODE_BIOMATRIX, _root.child(XML_NODE_RAWRES));
        for(auto& m: bioMatrix.matrices) {
            pugi::xml_document doc;
            doc.load_string(_biomatrix_xml.c_str());
            pugi::xml_node node = doc.child(XML_NODE_MATRIX);
            node.attribute(XML_ATTR_ID) = to_string(m.id).c_str();
            node.child(XML_NODE_APPLY_TO_CHAIN).attribute(XML_ATTR_CHAINID) = m.sourceChainId.c_str();
            node.child(XML_NODE_APPLY_TO_CHAIN).attribute(XML_ATTR_NEW_CHAINID) = m.newChainId.c_str();
            setTMatrix(node, m.tmatrix);
            biom_node.append_copy(node);
        }

    }

    vector<Tmdet::ValueObjects::Membrane> Xml::getMembranes() {
        vector<Tmdet::ValueObjects::Membrane> membranes;
        for (pugi::xml_node m_node = _root.child(XML_NODE_MEMBRANE); m_node; m_node = m_node.next_sibling(XML_NODE_MEMBRANE)) {
            pugi::xml_node tnode = m_node.child(XML_NODE_TMATRIX);
            Tmdet::ValueObjects::Membrane m = {
                getTMatrix(tnode),
                m_node.child(XML_NODE_NORMAL).attribute(XML_ATTR_Z).as_double(),
                0.0,
                0.0,
                Tmdet::Types::Membranes.at("Plain")
            };
            membranes.emplace_back(m);
        }
        return membranes;
    }

    void Xml::setMembranes(vector<Tmdet::ValueObjects::Membrane>& membranes) {
        
        for(auto& m : membranes) {
            pugi::xml_document doc;
            doc.load_string(_membrane_xml.c_str());
            pugi::xml_node node = doc.child(XML_NODE_MEMBRANE);
            node.child(XML_NODE_NORMAL).attribute(XML_ATTR_Z) = to_string(m.h).c_str();
            setTMatrix(node,m.tmatrix);
            _root.append_copy(node);
        }
    }

    void Xml::getChains(vector<Tmdet::ValueObjects::Chain>& chains) {
        for (pugi::xml_node c_node = _root.child(XML_NODE_CHAIN); c_node; c_node = c_node.next_sibling(XML_NODE_CHAIN)) {
            string type = c_node.attribute(XML_ATTR_TYPE).as_string();
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
                cerr << "Could not find gemmi chain for " << c_node.attribute(XML_ATTR_CHAINID).as_string() << endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    void Xml::setChains(vector<Tmdet::ValueObjects::Chain>& chains) {
        for(auto& c: chains) {
            pugi::xml_document doc;
            doc.load_string(_chain_xml.c_str());
            pugi::xml_node node = doc.child(XML_NODE_CHAIN);
            node.attribute(XML_ATTR_CHAINID) = c.id.c_str();
            node.attribute(XML_ATTR_NUM_TM) = to_string(c.numtm).c_str();
            node.attribute(XML_ATTR_TYPE) = c.type.name.c_str();
            node.child(XML_NODE_SEQ).text() = c.seq.c_str();
            setRegions(node,c.regions);
            _root.append_copy(node);
        }
    }

    vector<Tmdet::ValueObjects::Region> Xml::getRegions(pugi::xml_node& cnode) {
        vector<Tmdet::ValueObjects::Region> regions;
        for (pugi::xml_node r_node = cnode.child(XML_NODE_REGION); r_node; r_node = r_node.next_sibling(XML_NODE_REGION)) {
            char type = r_node.attribute(XML_ATTR_type).value()[0];
            Tmdet::ValueObjects::Region r = {
                r_node.attribute(XML_ATTR_SEQ_BEG).as_int(),
                r_node.attribute(XML_ATTR_PDB_BEG).as_string(),
                r_node.attribute(XML_ATTR_SEQ_END).as_int(),
                r_node.attribute(XML_ATTR_PDB_END).as_string(),
                0,
                0,
                Tmdet::Types::Regions.at(type)
            };
            regions.emplace_back(r);
        }
        return regions;
    }

    void Xml::setRegions(pugi::xml_node& pnode, vector<Tmdet::ValueObjects::Region>& regions) {
        for(auto& r: regions) {
            cerr << r.type.code << endl;
            char s[2];
            sprintf(s,"%c",r.type.code);
            pugi::xml_node node = pnode.append_child(XML_NODE_REGION);
            node.append_attribute(XML_ATTR_SEQ_BEG) = to_string(r.beg).c_str();
            node.append_attribute(XML_ATTR_PDB_BEG) = r.begi.c_str();
            node.append_attribute(XML_ATTR_SEQ_END) = to_string(r.end).c_str();
            node.append_attribute(XML_ATTR_PDB_END) = r.endi.c_str();
            node.append_attribute(XML_ATTR_type) = s;
        }
    }

}