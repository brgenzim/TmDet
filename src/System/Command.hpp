// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>

/**
 * @brief namespace for tmdet system
 *
 * @namespace Tmdet
 * @namespace System
 */
namespace Tmdet::System {

    /**
     * Executing external command
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
