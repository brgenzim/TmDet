#include <string>
#include <Version.hpp>
#include <Config.hpp>
#include <System/Date.hpp>
#include <Types/Protein.hpp>
#include <VOs/Xml.hpp>

namespace Tmdet::VOs {

    void Xml::notTransmembrane() {
        tmp = false;
        version = (version==""?Tmdet::version():version);
        modifications.emplace_back(
            Tmdet::System::Date::get(),
            (std::string)"Not transmembrane protein"
        );
        type = Tmdet::Types::ProteinType::SOLUBLE;
        bioMatrix.matrices.clear();
        bioMatrix.deletedChainIds.clear();
        membranes.clear();
        chains.clear();
    }

    
}
