#pragma once

#include <unordered_map>
#include <string>

namespace Tmdet::System {
    
    struct _arg {
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
            void define(bool mandatory, std::string shortFlag, std::string longFlag,
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
    };
    
}
