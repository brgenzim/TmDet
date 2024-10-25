#include <format>
#include <fstream>
#include <string>
#include <zlib.h>
#include <Exceptions/IOException.hpp>
#include <Helpers/Gzip.hpp>

namespace Tmdet::Helpers::Gzip {

    void compress(const std::string& source, const std::string& destination) {
        std::ifstream input(source, std::ios_base::binary | std::ios::ate);
        std::ofstream output(destination, std::ios_base::binary);

        auto size = input.tellg();
        input.seekg(0);
        std::string content(size, '\0');
        if (input.read(&content[0], size)) {
            writeFile(destination, content);
        }
    }

    void writeFile(const std::string& destination, const std::string& content) {
        gzFile file = gzopen(destination.c_str(), "wb9");
        if (file == nullptr) {
            throw Tmdet::Exceptions::IOException(std::format("Could not open '{}'", destination));
        }
        gzwrite(file, &content[0], content.length());
        int errorNum = 0;
        auto* message = gzerror(file, &errorNum);
        if (errorNum != Z_OK) {
            throw Tmdet::Exceptions::IOException(
                std::format("gzwrite failed: '{}', file name: '{}'", message, destination)
            );
        }

        gzclose(file);
        message = gzerror(file, &errorNum);
        if (errorNum != Z_OK) {
            throw Tmdet::Exceptions::IOException(
                std::format("gzclose failed: '{}', file name: '{}'", message, destination)
            );
        }
    }

}
