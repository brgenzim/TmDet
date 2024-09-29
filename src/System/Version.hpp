#ifndef __TMDET_SYSTEM_VERSION__
#define __TMDET_SYSTEM_VERSION__

#include <string>
#include <System/Command.hpp>

namespace Tmdet::System {

    /**
     * @brief Helper for getting version information
     * 
     */
    class Version {

        public:
            
            /**
             * @brief get the git version of the program as a string
             * 
             * @return std::string 
             */
            static std::string get() {
                return Command::run("git describe --tags --abbrev=0");
            }
    };
}

#endif