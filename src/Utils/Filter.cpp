// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt


#include <iostream>
#include <filesystem>
#include <vector>
#include <Helpers/String.hpp>
#include <System/Command.hpp>
#include <System/FilePaths.hpp>
#include <Utils/Filter.hpp>
#include <Utils/Md5.hpp>
#include <VOs/Chain.hpp>

namespace Tmdet::Utils {

    void Filter::run(bool applyTmFilter) {
        methodsDir = environment.get("CCTOP_METHODS_ROOT",DEFAULT_CCTOP_METHODS_ROOT);
        protein.tmFilterResults = false;
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                if (applyTmFilter && !chain.hasUnknownResidue()) {
                    if (createTempFasta(chain)) {
                        int tmp = 0;
                        runSignalP(chain.id);
                        parseSignalP(chain);
                        tmp += runPhobius(chain.id);
                        tmp += runScampi(chain.id);
                        tmp += runTMHMM(chain.id);
                        DEBUG_LOG("Filter results in {}",tmp);
                        chain.isTmp = (tmp>1);
                        protein.tmFilterResults |= chain.isTmp;
                    }
                }
                else {
                    chain.isTmp = true;
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
        bool ret = filePutContents(path,std::format(">{}_{}\n{}",protein.code,chain.id,chain.seq.substr(chain.signalP[1])));
        path = tempDir + "/list";
        ret &= filePutContents(path,std::format("{}/{}_{}.fas",tempDir,protein.code,chain.id));
        return true;
    }

    void Filter::runSignalP(std::string id) {
        std::string cmd = std::format("/usr/local/bin/signalp6 -fasta {}/{}_{}.fas -org euk --mode fast -od {}",
            tempDir,protein.code,id,tempDir);
        Tmdet::System::Command::run(cmd);
    }

    void Filter::parseSignalP(Tmdet::VOs::Chain& chain) {
        std::string path = tempDir + "/output.gff3";
        std::ifstream file(path);
        if (!file.is_open()) {
            logger.warn("Could not read file. Path: {}",path);
            return;
        }
        std::string line;
        std::getline(file, line);
        if (file.eof()) {
            return;
        }
        std::getline(file, line);
        
        auto fields = Tmdet::Helpers::String::explode("\t",line);
        if (fields.size() > 4 && fields[2] == "signal_peptide") {
            chain.signalP[0] = std::stoi(fields[3]);
            chain.signalP[1] = std::stoi(fields[4]);
            createTempFasta(chain);
            for (int i=chain.signalP[0]-1; i<chain.signalP[1]; i++) {
                chain.residues[i].selected = false;
            }
        }
        file.close();
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
