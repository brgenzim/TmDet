#include <fstream>
#include <iostream>
#include <vector>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include <gemmi/util.hpp>
#include <gemmi/math.hpp>

namespace cif = gemmi::cif;

__attribute__((noreturn)) void usage(const char *fileName);
gemmi::Vec3 selectAtoms(cif::Document doc, char *chid, int start, int end);


int main(int argc, char *argv[]) {
    int s1=-1;
    int e1=-1;
    int s2=-1;
    int e2=-1;
    auto input = (char *)nullptr;
    char chid = '-';
    for(int i=1; i<argc-1; i++) {
        if (!strcmp(argv[i],"-i")) {
            input = argv[i+1];
        }
        if (!strcmp(argv[i],"-s1")) {
            s1 = atoi(argv[i+1]);
        }
        if (!strcmp(argv[i],"-e1")) {
            e1 = atoi(argv[i+1]);
        }
        if (!strcmp(argv[i],"-s2")) {
            s2 = atoi(argv[i+1]);
        }
        if (!strcmp(argv[i],"-e2")) {
            e2 = atoi(argv[i+1]);
        }
    }
    if (input == (char *)nullptr || s1 == -1 || s2 == -1 || e1 == -1 || e2 == -1) {
        usage(argv[0]);
    }
    cif::Document doc = cif::read(gemmi::MaybeGzipped(input));
    auto center1 = selectAtoms(doc, &chid, s1, e1);
    auto center2 = selectAtoms(doc, &chid, s2, e2);
    auto dist = center1.dist(center2);
    std::cout << dist << std::endl;
}

void usage(const char *fileName) {
    std::cerr << "Usage: " << fileName << R"(
    -i : input cif file
)";
    exit(EXIT_FAILURE);
}

gemmi::Vec3 selectAtoms(cif::Document doc, char *chid, int start, int end) {
    auto mean = gemmi::Vec3(0.0,0.0,0.0);
    unsigned int count = 0;
    mean.x = mean.y = mean.z = 0;
    for (auto row : doc.sole_block().find("_atom_site.",
                          {"group_PDB", "label_asym_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z"})) {
        if (row[0] == "ATOM" &&
            atoi(row[2].c_str()) >= start &&
            atoi(row[2].c_str()) <= end) {
                if (*chid == '-') {
                    *chid = row[1].c_str()[0];
                }
                if (*chid == row[1].c_str()[0]) {
                    mean.x += atof(row[3].c_str());
                    mean.y += atof(row[4].c_str());
                    mean.z += atof(row[5].c_str());
                    //std::cout << row[2] << " :: ";
                    //std::cout << row[3] << " :: " << row[4] << " :: " << row[5] << " ==> ";
                    //std::cout << mean.x << " :: " << mean.y << " :: " << mean.z << std::endl;
                    count++;
                }
        }
    }
    mean /= count;
    return mean;
}

