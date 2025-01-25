// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt


#include <iostream>
#include <Utils/Filter.hpp>
#include <Utils/Md5.hpp>
#include <VOs/Chain.hpp>

namespace Tmdet::Utils {

    bool Filter::run() {
        bool ret = false;
        methodsDir = environment.get("CCTOP_METHODS_ROOT",DEFAULT_CCTOP_METHODS_ROOT);
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                if (createTempFasta(chain)) {
                    int tmp = 0;
                    tmp += runPhobius(chain.id);
                    tmp += runScampi(chain.id);
                    tmp += runTMHMM(chain.id);
                    chain.isTmp = (tmp>0);
                    ret |= chain.isTmp;
                }
            }
        );
        return ret;
    }

    bool Filter::filePutContents(std::string filePath, std::string content) {
        std::ofstream file(filePath);
        if (!file.is_open()) {
            logger.warn("Could not write file. Path: {}",filePath);
            return false;
        }
        file << content << std::endl;
        file.close();
        return true;
    }

    bool Filter::createTempFasta(Tmdet::VOs::Chain& chain ) {
        void* sig = hashing::md5::hash(raw);
        std::string hash =  hashing::md5::sig2hex(sig);
        tempDir = Tmdet::System::FilePaths::temp(hash);
        std::filesystem::create_directories(dir);
        std::string path = dir + "/" + protein.code + "_" + chain.id + ".fas";
        bool ret = filePutContents(path,std::format(">{}_{}\n{}",protein.code,chain.id,chain.seq));
        path = dir + "/list";
        ret &= filePutContents(path,std::format("{}_{}.fas",protein.code,chain.id));
        return true;
    }

    int Filter::runPhobius(std::string id) {
        std::string cmd = std::format("{}/phobius/phobius.pl {}/{}_{}",
                    methodsDir,tempDir,protein.code,id);
        return parsePhobius(Tmdet::System::Command::run(cmd));
    }

    int Filter::runScampi(std::string id) {
        std::string cmd = std::format("{}/bin/modhmms_scampi -f fa -s {}/list -m {}/scampi/DGHMM_KR_21.txt -o {} -r {}/scampi/replacement_letter_multi.rpl --nopostout --viterbi -u -L",
                    methodsDir,tempDir, methodsDir, tempDir, methodsDir);
        return parseScampi(Tmdet::System::Command::run(cmd));
    }

    int Filter::runTMHMM(std::string id) {
        std::string cmd = std::format("{}/tmhmm-2.0c/bin/decodeanhmm.Linux_x86_64 -modelfile {}/tmhmm-2.0c/lib/TMHMM2.0.model < {}/{}_{}.fas",
                    methodsDir, methodsDir, tempDir, protein.code, id);
        return parseTMHMM(Tmdet::System::Command::run(cmd));
    }

    int parsePhobius(std::string results) {
        int ret = 0;

        return ret;
    }

    int parseScampi(std::string results) {
        int ret = 0;

        return ret;
    }

    int parseTMHMM(std::string results) {
        int ret = 0;

        return ret;
    }
}
