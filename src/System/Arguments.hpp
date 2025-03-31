// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <unordered_map>
#include <string>

/**
 * @brief namespace for tmdet system
 *
 * @namespace Tmdet
 * @namespace System
 */
namespace Tmdet::System {
    
    /**
     * @brief container for arguments
     */
    struct _arg {
        /**
         * @brief flag to define if the argument is mandatory or not
         */
        bool mandatory;

        bool show;

        /**
         * @brief flag that gives if the argument is actually given in
         *        the command line
         */
        bool has;

        /**
         * @brief short form of the argument (like -h)
         */
        std::string shortFlag;

        /**
         * @brief long form of the argument (like --help)
         */
        std::string longFlag;

        /**
         * @brief description of the argument
         */
        std::string description;

        /**
         * @brief type of the argument (bool, string, int)
         */
        std::string type;

        /**
         * @brief default value of the argument (if it is not
         *        given in the command line)
         */
        std::string defaultValue;

        /**
         * @brief actual value of the argument (given in the
         *        command line)
         */
        std::string value;
    };

    /**
     * @brief handling command line arguments
     *
     */
    class Arguments {
        
        private:
            /**
             * @var list of arguments
             */
            std::unordered_map<std::string, _arg> _args;

            /**
             * @brief 
             * 
             * @param flag 
             * @return string 
             */
            std::string _setValueByLongFlag(char *flag);

            /**
             * @brief 
             * 
             * @param flag 
             * @return string 
             */
            std::string _setValueByShortFlag(char *flag);

            /**
             * @brief all arguments in a concatenated string
             */
            std::string commandLine = "";

            /**
             * @brief 
             * 
             * @param name 
             * @param value 
             */
            void _setValue(std::string name, char *value);

        public:

            /**
             * @brief Construct a new Arguments object
             * 
             */
            explicit Arguments();
            
            /**
             * @brief define argument
             *
             * @param mandatory : flag if the argument is mandatory
             * @param shortFlag : short flag, without '-'
             * @param longFlag : long flag, without '--'
             * @param descr : description of the argument
             * @param type : type, can be bool, int, real, double, string
             * @param defaultValue : default value (if not mandatory)
             */
            void define(bool mandatory, bool show, std::string shortFlag, std::string longFlag,
                std::string descr, std::string type, std::string defaultValue);

            /**
             * @brief get command line arguments
             *
             * @param argc
             * @param *argv[]
             * @return void
             */
            void set(int argc, char *argv[]);

            /**
             * @brief check arguments
             *
             * @return void or exit
             */
            void check();

            /**
             * @brief list arguments
             *
             * @return void
             */
            void list();
            
            /**
             * @brief get value as bool of an argument
             *
             * @param name
             * @return bool
             */
            bool getValueAsBool(std::string name);

            /**
             * @brief get value as int of an argument
             *
             * @param name
             * @return int
             */
            int getValueAsInt(std::string name);

            /**
             * @brief get value as string of an argument
             *
             * @param name
             * @return std::string
             */
            std::string getValueAsString(std::string name);

            /**
             * @brief get value as float of an argument
             *
             * @param name
             * @return float
             */
            float getValueAsFloat(std::string name);

            /**
             * @brief return the concatenated argument string
             * 
             * @return std::string 
             */
            std::string getCommandLine() const;

            void setCommandLine();

            std::string hidePaths(std::string input);
    };
    
}
