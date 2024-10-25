#pragma once

#include<string>

namespace Tmdet::Helpers::Gzip {

    /**
     * @brief compress content of source into destination.
     * @attention This function is created for small CIF files intentionally;
     *            source file must be small enough for available memory.
     */
    void compress(const std::string& source, const std::string& destination);

    /**
     * @brief Gzip and write data into destination file.
     */
    void writeFile(const std::string& destination, const std::string& data);

}
