#pragma once

#include <string>
#include <ios>

namespace Tmdet::Exceptions {

    class IOException : public std::ios_base::failure {
    public:
        IOException(const std::string& message) : std::ios_base::failure(message) {}
        IOException(const char* message) : std::ios_base::failure(message) {}
    };

}
