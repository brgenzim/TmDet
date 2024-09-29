#ifndef __TMDET_SYSTEM_COMMAND__
#define __TMDET_SYSTEM_COMMAND__

#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>

namespace Tmdet::System {

    /**
     * Helper for executing external command
     */
    class Command {

        public:

            /**
             * @brief run other program and get std output as string
             * 
             * @param cmd 
             * @return std::string 
             */
            static std::string run( const std::string &cmd ) {
                FILE *pfd = popen( cmd.c_str(), "r" );
                std::string res="";
                if ( pfd != nullptr ) {    
                    while ( !feof(pfd) ) {
                        char buf[ 1024 ] = {0};
                        if ( fgets(buf, sizeof(buf), pfd) != (char *)nullptr ) {
                            res += ( std::string(buf) );
                        }
                    }
                    pclose( pfd );
                }
                return res;
            }
    };
}

#endif
