// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <pugixml.hpp>
#include <DTOs/XmlRW/Constants4.hpp>
#include <DTOs/XmlRW/Reader.hpp>
#include <Exceptions/SyntaxErrorException.hpp>

namespace Tmdet::DTOs::XmlRW {

    void Reader::read(const std::string& path) {
        if (pugi::xml_parse_result result = doc.load_file(path.c_str()); !result) {
            throw Tmdet::Exceptions::SyntaxErrorException(path,(int)result.offset,result.description());
        }
    }

    bool Reader::isVersion4() {
        auto root = doc.child(XML_NODE_ROOT);
        return (root.child(XML_NODE_RAWDATA)?true:false);
    }

    void Reader::readXml(Tmdet::VOs::Xml& xmlData, const std::string& path) {
        read(path);
        if (isVersion4()) {
            reader4.setRoot(doc);
            reader4.readXml(xmlData);
        }
        else {
            reader3.setRoot(doc);
            reader3.readXml(xmlData);
        }
    }
}