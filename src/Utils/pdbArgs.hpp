#ifndef __UNITMP_PDBLIB_UTILS_PDBARGS__
#define __UNITMP_PDBLIB_UTILS_PDBARGS__

#include <unordered_map>
#include <string>

namespace UniTmp::PdbLib::Utils {
    
    struct _pdbArg {
        bool mandatory;
        bool has;
        std::string shortFlag;
        std::string longFlag;
        std::string description;
        std::string type;
        std::string defaultValue;
        std::string value;
    };

    /**
     * handling command line arguments
     */
    class pdbArgs {
        
        private:
            std::unordered_map<std::string, _pdbArg> _args;

            std::string _setValueByLongFlag(char *flag);
            std::string _setValueByShortFlag(char *flag);
            void _setValue(std::string name, char *value);

        public:
            pdbArgs();
            
            /**
             * define argument
             * @param bool mandatory : flag if the argument is mandatory
             * @param std::string *shortFlag : short flag, without '-'
             * @param std::string longFlag : long flag, without '--'
             * @param std::string descr : description of the argument
             * @param std::string type : type, can be bool, int, real, double, string
             * @param std::string defaultValue : default value (if not mandatory)
             * @return void
             */
            void define(bool mandatory, std::string shortFlag, std::string longFlag,
                std::string descr, std::string type, std::string defaultValue);

            /**
             * get command line arguments
             * @param int argc
             * @param char *argv[]
             * @return void
             */
            void set(int argc, char *argv[]);

            /**
             * check arguments
             * @return void or exit
             */
            void check();

            /**
             * list arguments
             * @return void
             */
            void list();
            
            /**
             * get value as bool of an argument
             * @param std::string name
             * @return bool
             */
            bool getValueAsBool(std::string name);

            /**
             * get value as int of an argument
             * @param std::string name
             * @return int
             */
            int getValueAsInt(std::string name);

            /**
             * get value as string of an argument
             * @param std::string name
             * @return std::string
             */
            std::string getValueAsString(std::string name);
    };
    
}
#endif