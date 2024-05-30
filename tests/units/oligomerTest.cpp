#include <iostream>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>
#include <Utils/Oligomer.hpp>

#define PDB_CIF_PATH "/home/data/wwPDB/data/structures/divided/updated_mmcif/"
#define IS_TRUE(x, n) { \
        if ((x)) { \
            cerr << "Passed: " << __FUNCTION__ << endl; \
        } \
        else { \
            cout << __FUNCTION__ << " failed on line " << __LINE__ << " type: " << n << endl; \
        } \
}


using namespace std;

Tmdet::Types::Oligomer getType(string code) {
    string path = PDB_CIF_PATH + code.substr(1,2) + "/" + code + "_updated.cif.gz";
    gemmi::cif::Document doc = gemmi::cif::read(gemmi::MaybeGzipped(path));
    return Tmdet::Utils::Oligomer::getOligomerType(doc);
}

void checkMonomer() {
    string type = getType("2guf").name;
    IS_TRUE( (type == Tmdet::Types::OligomerType::MONOMER.name), type);
}

void checkHomoOligomer() {
    string type = getType("1a0s").name;
    IS_TRUE( (type == Tmdet::Types::OligomerType::HOMO_OLIGOMER.name), type);
}

void checkHomoHeteroOligomer() {
    string type = getType("1occ").name;
    IS_TRUE( (type == Tmdet::Types::OligomerType::HOMO_HETERO_OLIGOMER.name), type);
}

void checkHeteroWithHomoOligomer() {
    string type = getType("8j0s").name;
    IS_TRUE( (type == Tmdet::Types::OligomerType::HETERO_WITH_HOMO_OLIGOMER.name), type);
}

void checkHeteroOligomer() {
    string type = getType("8jpn").name;
    IS_TRUE( (type == Tmdet::Types::OligomerType::HETERO_OLIGOMER.name), type);
}

int main(void) {
    checkMonomer();
    checkHomoOligomer();
    checkHomoHeteroOligomer();
    checkHeteroWithHomoOligomer();
    checkHeteroOligomer();
}