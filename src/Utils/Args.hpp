#ifndef __TMDET_UTILS_ARGS__
#define __TMDET_UTILS_ARGS__

#include <unordered_map>
#include <string>

using namespace std;

namespace Tmdet::Utils {
    
    struct _arg {
        bool mandatory;
        bool has;
        string shortFlag;
        string longFlag;
        string description;
        string type;
        string defaultValue;
        string value;
    };

    /**
     * handling command line arguments
     */
    class Args {
        
        private:
            unordered_map<string, _arg> _args;

            string _setValueByLongFlag(char *flag);
            string _setValueByShortFlag(char *flag);
            void _setValue(string name, char *value);

        public:
            Args();
            
            /**
             * define argument
             * @param bool mandatory : flag if the argument is mandatory
             * @param string *shortFlag : short flag, without '-'
             * @param string longFlag : long flag, without '--'
             * @param string descr : description of the argument
             * @param string type : type, can be bool, int, real, double, string
             * @param string defaultValue : default value (if not mandatory)
             * @return void
             */
            void define(bool mandatory, string shortFlag, string longFlag,
                string descr, string type, string defaultValue);

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
             * @param string name
             * @return bool
             */
            bool getValueAsBool(string name);

            /**
             * get value as int of an argument
             * @param string name
             * @return int
             */
            int getValueAsInt(string name);

            /**
             * get value as string of an argument
             * @param string name
             * @return string
             */
            string getValueAsString(string name);
    };
    
}
#endif