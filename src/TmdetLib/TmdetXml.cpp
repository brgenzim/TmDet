#include <string>
#include <vector>
#include <iostream>
#include <pugixml.hpp>
#include <TmdetXml.hpp>
#include <TmdetStruct.hpp>

using namespace std;
using namespace gemmi;

namespace UniTmp::TmdetLib {

    TmdetXml::TmdetXml() {

    }

    TmdetXml::~TmdetXml() {

    }
    
    void TmdetXml::read(string path) {
        pugi::xml_parse_result result = _doc.load_file(path.c_str());
        if (!result) {
            cerr << "Error description: " << result.description() << endl;
            cerr << "Error offset: " << result.offset << endl;
            exit(EXIT_FAILURE);
        }
        _root = _doc.child(XML_NODE_ROOT);
    }

    void TmdetXml::write(string path) {
        _doc.save_file(path.c_str(),"  ");
    }

    void TmdetXml::create() {
        _doc.load_string(_pdbtm_xml.c_str());
        _root = _doc.child(XML_NODE_ROOT);
    }

    bool TmdetXml::getTmp() {
        return _root.attribute(XML_ATTR_TMP).as_bool();
    }

    void TmdetXml::setTmp(bool tmp) {
        _root.attribute(XML_ATTR_TMP).set_value(tmp?"yes":"no");
    }

    string TmdetXml::getCode() {
        return _root.attribute(XML_ATTR_ID).value();
    }

    void TmdetXml::setCode(string code) {
        _root.attribute(XML_ATTR_ID).set_value(code.c_str());
    }

    string TmdetXml::getCreateDate() {
        return _root.child(XML_NODE_CREATE_DATE).text().get();
    }

    void TmdetXml::setCreateDate(string date) {
        _root.child(XML_NODE_CREATE_DATE).text() = date.c_str();
    }

    vector<_tmdetModification> TmdetXml::getModifications() {
        vector<_tmdetModification> mods;
        for (pugi::xml_node mod = _root.child(XML_NODE_MODIFICATION); mod; mod = mod.next_sibling(XML_NODE_MODIFICATION)) {
            _tmdetModification m = {mod.child(XML_NODE_DATE).text().get(),mod.child(XML_NODE_DESCR).text().get()};
            mods.emplace_back(m);
        }
        return mods;
    }

    void TmdetXml::setModifications(vector<_tmdetModification> mods) {
        pugi::xml_node node;
        for(auto m : mods) {
            node = _root.insert_child_after(XML_NODE_MODIFICATION, _root.child(XML_NODE_CREATE_DATE));
            pugi::xml_node date_node = node.append_child(XML_NODE_DATE);
            date_node.append_child(pugi::node_pcdata).set_value(m.date.c_str());
            pugi::xml_node descr_node = node.append_child(XML_NODE_DESCR);
            descr_node.append_child(pugi::node_pcdata).set_value(m.descr.c_str());
        }
    }

    double TmdetXml::getQvalue() {
        return _root.child(XML_NODE_RAWRES).child(XML_NODE_TMRES).text().as_double();
    }

    void TmdetXml::setQvalue(double q) {
        char s[100];
        sprintf(s,"%.2f",q);
        _root.child(XML_NODE_RAWRES).child(XML_NODE_TMRES).text() = s;
    }

    string TmdetXml::getTmtype() {
        return _root.child(XML_NODE_RAWRES).child(XML_NODE_TMTYPE).text().get();
    }

    void TmdetXml::setTmtype(const char* ptype) {
        cerr << "::" << ptype << "::" << endl;
        _root.child(XML_NODE_RAWRES).child(XML_NODE_TMTYPE).text() = ptype;
    }
}