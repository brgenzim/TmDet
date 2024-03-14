
#ifndef __UNITMP_PDBLIB_TYPES_PDB_MOLECULE_TYPES__
#define __UNITMP_PDBLIB_TYPES_PDB_MOLECULE_TYPES__

#include <unordered_map>

#include "PdbAtomTypes.hpp"
using namespace UniTmp::PdbLib::Types::PdbAtomTypes;

namespace UniTmp::PdbLib::Types {

    struct PdbMoleculeType {
        const char *name;
        const char *a1;
        const char *a3;
        std::unordered_map<const char*, PdbAtomType> atoms;
    };

    namespace PdbMoleculeTypes {
        const PdbMoleculeType ALA {
                "Alanine", "ALA", "A", 
                {{"N", N},{"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"OXT", O}}
        }; 
        const PdbMoleculeType CYS {
                "Cysteine", "CYS", "C",
                {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"SG", S},{"OXT", O}}
        };
        const PdbMoleculeType ASP {
                "Aspartic acid", "ASP", "D",
                {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_CAR}, {"OD1", O}, {"OD2", O}, {"OXT", O}} 
        };
        const PdbMoleculeType GLU {
                "Glutamic acid", "GLU", "E",
                {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_ALI}, {"CD", C_CAR}, {"OE1", O}, {"OE2", O}, {"OXT", O}} 
        };
        const PdbMoleculeType PHE {
            "Phenylalanine", "PHE", "F",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_CAR}, {"CD1", C_CAR}, {"CD2", C_CAR}, {"CE1", C_CAR}, {"CE2", C_CAR}, {"CZ", C_CAR}, {"OXT", O}} 
        };
        const PdbMoleculeType GLY {
            "Glycine", "GLY", "G", 
            {{"N", N},{"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"OXT", O}}
        };
        const PdbMoleculeType HIS {
            "Histidine", "HIS", "H",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_CAR}, {"ND1", N_AMD}, {"CD2", C_CAR}, {"CE1", C_CAR}, {"NE2", N_AMD}, {"OXT", O}} 
        };
        const PdbMoleculeType ILE {
            "Isoleucine", "ILE", "I",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG1", C_ALI}, {"CG2", C_ALI}, {"CD1", C_ALI}, {"OXT", O}} 
        };
        const PdbMoleculeType LYS {
            "Lysine", "LYS", "K",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_ALI}, {"CD", C_ALI}, {"CE", C_ALI}, {"NZ", N_AMN}, {"OXT", O}} 
        };
        const PdbMoleculeType LEU {
            "Leucine", "LEU", "L",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_ALI}, {"CD1", C_ALI}, {"CD2", C_ALI}, {"OXT", O}}
        };
        const PdbMoleculeType MET {
            "Methionine", "MET", "M",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_ALI}, {"SD", S}, {"CE", C_ALI}, {"OXT", O}}
        };
        const PdbMoleculeType ASN {
            "Asparagine", "ASN", "N",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_CAR}, {"OD1", O}, {"ND2", N_AMD}, {"OXT", O}}
        };
        const PdbMoleculeType PRO {
            "Proline", "PRO", "P",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_ALI}, {"CD", C_ALI}, {"OXT", O}}
        };
        const PdbMoleculeType GLN {
            "Glutamine", "GLN", "Q",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_ALI}, {"CD", C_CAR}, {"OE1", O}, {"NE2", N_AMD}, {"OXT", O}}
        };
        const PdbMoleculeType ARG {
            "Arginine", "ARG", "R",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_ALI}, {"CD", C_ALI}, {"NE", N_AMD}, {"CZ", C_CAR}, {"NH1", N_AMD}, {"NH2", N_AMD}, {"OXT", O}}
        };
        const PdbMoleculeType SER {
            "Serine", "SER", "S",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_ALI}, {"OG", O}, {"OXT", O}}
        };
        const PdbMoleculeType THR {
            "Threonine", "THR", "T",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_ALI}, {"OG1", O}, {"CG2", C_ALI}, {"OXT", O}}
        };
        const PdbMoleculeType VAL {
            "Valine", "VAL", "V",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG1", C_ALI}, {"CG2", C_ALI}, {"OXT", O}}
        };
        const PdbMoleculeType TRP {
            "Tryptophan", "TRP", "W",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_CAR}, {"CD1", C_CAR}, {"CD2", C_CAR}, {"NE1", N_AMD}, {"CE2", C_CAR}, {"CE3", C_CAR}, {"CZ2", C_CAR}, {"CZ3", C_CAR}, {"CH2", C_CAR}, {"OXT", O}}
        };
        const PdbMoleculeType TYR {
            "Tyrosine", "TYR", "Y",
            {{"N", N}, {"CA", C_ALI}, {"C", C_CAR}, {"O", O}, {"CB", C_ALI}, {"CG", C_CAR}, {"CD1", C_CAR}, {"CD2", C_CAR}, {"CE1", C_CAR}, {"CE2", C_CAR}, {"CZ", C_CAR}, {"OH", O}, {"OXT", O}}
        };
        const PdbMoleculeType ACE {
            "Acetyl group", "ACE", "",
            {{"C", C}, {"O", O}, {"CH3", C_ALI}}
        };
        const PdbMoleculeType NME {
            "Methylamine", "NME", "",
            {{"N", N}, {"C", C}}
        };
        const PdbMoleculeType FOR {
            "Formyl group", "FOR", "",
            {{"C", C}, {"O", O}}
        };
        const PdbMoleculeType HOH {
            "Water", "HOH", "",
            {{"O", O}}
        };

        const std::unordered_map<const char*, PdbMoleculeType> stdAminoAcids {
            {"ALA", ALA}, {"CYS", CYS}, {"ASP", ASP}, {"GLU", GLU}, {"PHE", PHE},
            {"GLY", GLY}, {"HIS", HIS}, {"ILE", ILE}, {"LYS", LYS}, {"LEU", LEU},
            {"MET", MET}, {"ASN", ASN}, {"PRO", PRO}, {"GLN", GLN}, {"ARG", ARG},
            {"SER", SER}, {"THR", THR}, {"VAL", VAL}, {"TRP", TRP}, {"TYR", TYR}
        };

    }
}

#endif
