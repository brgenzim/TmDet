#include <string>
#include <pugixml.hpp>
#include <DTOs/Xml/BaseWriter.hpp>

namespace Tmdet::DTOs::Xml {

    void BaseWriter::write(const std::string& path) const {
        _doc.save_file(path.c_str(),"  ");
    }
}
