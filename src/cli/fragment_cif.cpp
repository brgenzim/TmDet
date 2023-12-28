#include <fstream>
#include <iostream>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/to_cif.hpp>

namespace cif = gemmi::cif;

void usage(char *fileName);
std::string createDir(const char *destDir, const char *cifName);
void writeCif(const char *cifPath, const cif::Block& block) ;

int main(int argc, char *argv[]) {

    char *input = (char *)NULL;
    char *destDir = (char *)NULL;
    bool subDir = false;
    for(int i=1; i<argc; i++) {
        if (!strcmp(argv[i],"-i")) {
            input = argv[++i];
        }
        if (!strcmp(argv[i],"-d")) {
            destDir = argv[++i];
        }
        if (!strcmp(argv[i],"-s")) {
            subDir = true;
        }
    }
    if (input == (char *)NULL || destDir == (char *)NULL) {
        usage(argv[0]);
    }
    cif::Document doc = cif::read(gemmi::MaybeGzipped(input));
    for (cif::Block& block : doc.blocks) {
        std::string cifPath = (subDir?
            createDir(destDir,block.name.c_str()):
            (std::string)destDir) + "/" + block.name + ".cif";
        writeCif(cifPath.c_str(), block);
        std::cerr << "Writing " << block.name << " to " << cifPath << std::endl;
    }
}

void usage(char *fileName) {
    std::cerr << R"(Usage:
    -i : input cif file
    -d : destination directory
    -s : create subdirectories
)";
    exit(EXIT_FAILURE);
}

std::string createDir(const char *destDir, const char *cifName) {
    std::string path = (std::string)destDir 
            + "/" + cifName[0]
            + "/" + cifName[1];
    std::string cmd = (std::string)"mkdir -p " + path;
    system(cmd.c_str());
    return path;
}

void writeCif(const char *cifPath, const cif::Block& block) {
    std::ofstream os;
    cif::WriteOptions wo;
    os.open(cifPath);
    write_cif_block_to_stream(os, block, wo);
    os.close();
}