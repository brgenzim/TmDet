// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt


#include <iostream>
#include <filesystem>
#include <System/Command.hpp>
#include <System/FilePaths.hpp>
#include <Utils/Filter.hpp>
#include <Utils/Md5.hpp>
#include <VOs/Chain.hpp>

namespace Tmdet::Utils {

    void Filter::run() {
        methodsDir = environment.get("CCTOP_METHODS_ROOT",DEFAULT_CCTOP_METHODS_ROOT);
        protein.tmFilterResults = false;
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                if (createTempFasta(chain)) {
                    int tmp = 0;
                    tmp += runPhobius(chain.id);
                    tmp += runScampi(chain.id);
                    tmp += runTMHMM(chain.id);
                    DEBUG_LOG("Filter results in {}",tmp);
                    chain.isTmp = (tmp>0);
                    protein.tmFilterResults |= chain.isTmp;
                }
            }
        );
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
        std::string hash =  Tmdet::Utils::Md5::getHash(chain.seq);
        tempDir = Tmdet::System::FilePaths::temp(hash);
        std::filesystem::create_directories(tempDir);
        std::string path = tempDir + "/" + protein.code + "_" + chain.id + ".fas";
        bool ret = filePutContents(path,std::format(">{}_{}\n{}",protein.code,chain.id,chain.seq));
        path = tempDir + "/list";
        ret &= filePutContents(path,std::format("{}/{}_{}.fas",tempDir,protein.code,chain.id));
        return true;
    }

    int Filter::runPhobius(std::string id) {
        std::string cmd = std::format("{}/phobius/phobius.pl {}/{}_{}.fas",
                    methodsDir,tempDir,protein.code,id);
        return parsePhobius(Tmdet::System::Command::run(cmd));
    }

    int Filter::runScampi(std::string id) {
        std::string cmd = std::format("{}/bin/modhmms_scampi -f fa -s {}/list -m {}/scampi/DGHMM_KR_21.txt -o {} -r {}/scampi/replacement_letter_multi.rpl --nopostout --viterbi -u -L",
                    methodsDir,tempDir, methodsDir, tempDir, methodsDir);
        return parseScampi(Tmdet::System::Command::run(cmd));
    }

    int Filter::runTMHMM(std::string id) {
        std::string cmd = std::format("{}/tmhmm-2.0c/bin/decodeanhmm.Linux_x86_64 -f {}/tmhmm-2.0c/lib/TMHMM2.0.options -modelfile {}/tmhmm-2.0c/lib/TMHMM2.0.model < {}/{}_{}.fas",
                    methodsDir, methodsDir,  methodsDir, tempDir, protein.code, id);
        return parseTMHMM(Tmdet::System::Command::run(cmd));
    }

    int Filter::parsePhobius(std::string results) {
        return (results.find("FT   TRANSMEM") == std::string::npos?0:1);
    }

    int Filter::parseScampi(std::string results) {
        return (results.find("<label>M</label>") == std::string::npos?0:1);
    }

    int Filter::parseTMHMM(std::string results) {
        int ret = 0;
        return (results.find("MMMMMMMMM") == std::string::npos?0:1);
        return ret;
    }
}
