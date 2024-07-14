#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>
#include <list>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <Utils/Alignment.hpp>
#include <DTOs/TmdetStruct.hpp>

using namespace std;

namespace Tmdet::Utils::Alignment {

    static gemmi::Residue* createResidue(int seqNum, int labelSeqNum, string name, string chainName);
    static bool compareResidues(const gemmi::Residue& res1, const gemmi::Residue& res2);
    static void initGaps(gemmi::Chain &chain, vector<int> &gaps);

#ifdef __ALIGNMENT_DBG
    void printMatrix(Eigen::MatrixXi& matrix);

#define MATRIX_DIR "/tmp/pdb_map"

struct alignState {

    Eigen::MatrixXi* scores;
    Eigen::MatrixXi* UTI;
    Eigen::MatrixXi* UTJ;
    int seqResNum; // num of residues from SEQRES/entity lines
    int atomResNum; // num of residues from atom lines
    int stateId; // a counter to help unique file name generation
    char fileName[200];
    string code;
    char chain;
};

void writeMatrix(FILE* file, Eigen::MatrixXi *mx, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        // write data
        for (int j = 0; j < cols; j++) {
            fprintf(file, "%5d", (*mx)(i, j));
            if (j+1 < cols) { fputc(' ', file); }
        }
        fputs("\n", file);
    }
}

void writeDebugState(alignState* state) {

    if (state->chain != 'C') {
        return;
    }


    char path[400];

    sprintf(state->fileName, "%s-align-matricies-%05d.txt",
        state->code.c_str(), state->stateId++);
    sprintf(path, "%s/%s", MATRIX_DIR, state->fileName);


    FILE* output = fopen(path, "w");
    if (!output) {
        perror("File opening failed");
        fprintf(stderr, "File name: %s\n", path);
        exit(1);
    }

    fprintf(output, "%s\n", path);
    fprintf(output, "scores:\n");
    writeMatrix(output, state->scores, state->seqResNum, state->atomResNum);

    fprintf(output, "UTI:\n");
    writeMatrix(output, state->UTI, state->seqResNum, state->atomResNum);
    fprintf(output, "UTJ:\n");
    writeMatrix(output, state->UTI, state->seqResNum, state->atomResNum);

    fclose(output);
}

#endif

    void alignSequences(gemmi::Chain &chain, vector<string> sequence) {

        vector<gemmi::Residue*> seqResidues;
        int seqNum = 1;
        for (auto &residueName : sequence) {
            auto residue = createResidue(seqNum, seqNum, residueName, chain.name);
            seqResidues.emplace_back(residue);
            seqNum++;
        }

        int alignmentStripe = 200; // formerly "alisav"
        int maxgap = 300;
        // number of residues based on atom lines:
        int chainLength = chain.residues.size();
        // number of residues based on SEQRES/entity_poly info:
        int seqResiduesLength = seqResidues.size();

        Eigen::MatrixXi scores(seqResiduesLength, chainLength);
        Eigen::MatrixXi iPath(seqResiduesLength, chainLength);
        Eigen::MatrixXi jPath(seqResiduesLength, chainLength);
        scores.setZero();
        iPath.setZero();
        jPath.setZero();

        if (abs(chainLength - seqResiduesLength) > 200) {
            alignmentStripe = abs(chainLength - seqResiduesLength) + 200;
        }
        if (abs(chainLength - seqResiduesLength) > 100) {
            maxgap = 1000;
        }
        vector<int> gaps;
        gaps.resize(chainLength);
        // gaps init
        initGaps(chain, gaps);

#ifdef __ALIGNMENT_DBG
        string code("5eit");
        alignState state;
        state.code = code;
        state.seqResNum = seqResiduesLength;
        state.atomResNum = chainLength;
        state.scores = &scores;
        state.UTI = &iPath;
        state.UTJ = &jPath;
        state.chain = chain.name[0];
        state.stateId = 0;
        writeDebugState(&state);
#endif


        // Calculating scores, iPath and jPath
        const int GAPOPEN = 2;
        int maxScore, i, j;
        for (i = 0; i < seqResiduesLength; i++) {
            for (j = (i-alignmentStripe < 0 ? 0 : i-alignmentStripe);
                (j < chainLength && j < i+alignmentStripe); j++) {

// #ifdef __ALIGNMENT_DBG
//                 if (chain.name[0] == 'C') {
//                     cout << "i: " << i << " j: " << j << endl;
//                 }
// #endif

                iPath(i, j) = i-1;
                jPath(i, j) = j-1;
                if (i == 0) { iPath(0, j) = -1; jPath(0, j) = j-1; }
                if (j == 0) { iPath(i, 0) = i-1; jPath(i, 0) = -1; }
                maxScore = 0;
                if (i > 0 && j > 0) { maxScore = scores(i-1, j-1); }
                if (i > 0) {
                    for (int k = j-2; (k >= 0 && k > j-maxgap); k--) {
                        int gap = 2 * (j-k) + GAPOPEN;
                        if (scores(i-1, k) - gap >= maxScore) {
                            maxScore = scores(i-1, k) - gap;
                            jPath(i, j) = k;
                        }
                    }
                }
                if (j > 0) {
                    for (int k = i-2; (k >= 0 && k > i-maxgap); k--) {
                        int gap = gaps[j-1] * (i-k) + GAPOPEN;
                        if (scores(k, j-1) - gap >= maxScore) {
                            maxScore = scores(k, j-1) - gap;
                            iPath(i, j) = k;
                        }
                    }
                }
                scores(i, j) = (seqResidues[i]->name == chain.residues[j].name ? 10 : 0) + maxScore;
#ifdef __ALIGNMENT_DBG
                writeDebugState(&state);
#endif
            }
        }
// #ifdef __ALIGNMENT_DBG
//         if (chain.name == "C") {
//             cout << "scores:" << endl; printMatrix(scores);
//             cout << "iPath:" << endl; printMatrix(iPath);
//             cout << "jPath:" << endl; printMatrix(jPath);
//         }
// #endif

        maxScore = scores(seqResiduesLength-1, chainLength-1);
        list<gemmi::Residue*> seqResiduesList(seqResidues.begin(), seqResidues.end());
        list<gemmi::Residue*> chainList;
        for (auto it = chain.residues.begin(); it < chain.residues.end(); it++) {
            chainList.push_back(&*it);
        }

        // trace back paths?
        int nuti = seqResiduesLength - 1;
        for (int i = nuti; i >= 0; i--) {
            if (maxScore <= scores(i, chainLength - 1)) {
                maxScore = scores(i, chainLength - 1);
                nuti = i;
            }
        }
        int nutj = chainLength - 1;
        for (int j = nutj; j >= 0; j--) {
            if (maxScore <= scores(seqResiduesLength - 1, j)) {
                maxScore = scores(seqResiduesLength - 1, j);
                nutj = j;
            }
        }

        // update seqResidues sequence
        for (int j = nutj+1; j < chainLength; j++) {
            auto residueFromAtomLine = chain.residues[j].empty_copy();
            residueFromAtomLine.atoms = chain.residues[j].atoms;
            seqResiduesList.push_back(&residueFromAtomLine);
        }

        // inserts
        int UTI = nuti;
        int UTJ = nutj;

        while (UTI >= 0 && UTJ >= 0) {
            auto currentResidue = chain.residues[UTJ];
            if (seqResidues[UTI]->name != currentResidue.name) {
                cerr << "Conflict in residue name: chain "
                    << chain.name << " at position " << currentResidue.label_seq.value
                    << seqResidues[UTI]->name << " vs " << currentResidue.name << endl;
                seqResidues[UTI]->name = currentResidue.name;
            }
            nuti = iPath(UTI, UTJ);
            nutj = jPath(UTI, UTJ);
            seqResidues[UTI]->seqid = chain.residues[UTJ].seqid;
            seqResidues[UTI]->label_seq = chain.residues[UTJ].label_seq;
            for (int j = UTJ-1; j > nutj; j--) {
                auto residueFromAtomLine = chain.residues[j].empty_copy();
                residueFromAtomLine.atoms = chain.residues[j].atoms;
                auto equals = [&seqResidues, &UTI](gemmi::Residue *item) { return item == seqResidues[UTI]; };
                auto pos = find_if(seqResiduesList.begin(), seqResiduesList.end(), equals);
                // insert pointer before position
                seqResiduesList.insert(pos, &residueFromAtomLine);
            }
            UTI = nuti;
            UTJ = nutj;
        }

#ifdef __ALIGNMENT_DBG
        printf("Chain %s: UTI: %d, UTJ: %d, nuti: %d, nutj: %d\n",
            chain.name.c_str(), UTI, UTJ, nuti, nutj);
#endif


        // std::vector<gemmi::Residue>  newChainResidues;
        // int labelSeqNum = 1;
        // for (auto& oneLetterSeq : withGaps) {

        //     if (oneLetterSeq != '-') {
        //         labelSeqNum++;
        //         continue;
        //     }

        //     gemmi::Residue* newResidue = createResidue(0, labelSeqNum, sequence[labelSeqNum - 1], chain.name);
        //     newChainResidues.emplace_back(*newResidue);
        //     labelSeqNum++;
        // }
        // chain.append_residues(newChainResidues);
        std::sort(chain.residues.begin(), chain.residues.end(), compareResidues);
    }

    gemmi::Residue* createResidue(int seqNum, int labelSeqNum, string name, string chainName) {
        auto residue = new gemmi::Residue();
        residue->seqid.num.value = seqNum;
        residue->label_seq.value = labelSeqNum;
        residue->name = name;
        residue->subchain = chainName;
        return residue;
    }

    void alignResidues(const Tmdet::ValueObjects::TmdetStruct& tmdetVO) {

        for(auto& chain: tmdetVO.gemmi.models[0].chains) {

            auto sequence = Tmdet::DTOS::TmdetStruct::getChainSequence(tmdetVO, chain);

            if (sequence.empty() || chain.residues.empty()) {
                // no supporting information to do gap fix
                // or there is no residues for the iteration
                continue;
            }

            alignSequences(chain, sequence);
        }
    }

    void initGaps(gemmi::Chain &chain, vector<int> &gaps) {
        const int GAP = 5;
        gaps[0] = 0;

        auto &residues = chain.residues;
        auto previous = residues.begin();
        auto current = previous + 1;
        for (int index = 1; current < residues.end(); current++, previous++, index++) {
            double atomDistance = 0.0;
            auto CA = previous->get_ca();
            auto currentCA = current->get_ca();
            if (CA != NULL && currentCA != NULL) {
                atomDistance = CA->pos.dist(currentCA->pos);
            }
            int seqDistance = current->label_seq.value - previous->label_seq.value;
            if (seqDistance > 1 && seqDistance < 500) {
                atomDistance = 20;
            }
            if (atomDistance > 18.3) {
                if (index <= GAP) {
                    gaps[index-1] = 2;
                } else {
                    gaps[index-1] = 0;
                }
                if (index >= (int)gaps.size() - GAP) {
                    gaps[index-1] = 2;
                }
            } else {
                gaps[index-1] = 2;
            }
        }
    }

    bool compareResidues(const gemmi::Residue& res1, const gemmi::Residue& res2) {
        return res1.label_seq.value < res2.label_seq.value;
    }

#ifdef __ALIGNMENT_DBG
    template <typename Iterator>
    void _printResiduesImpl(const Iterator &start, const Iterator &end);

    void printResidues(const vector<gemmi::Residue*> &residues) {
        _printResiduesImpl(residues.begin(), residues.end());
    }

    void printResidues(const list<gemmi::Residue*> &residues) {
        _printResiduesImpl(residues.begin(), residues.end());
    }

    template <typename Iterator>
    void _printResiduesImpl(const Iterator &start, const Iterator &end) {
        for_each(start, end,
            [](gemmi::Residue* item) { cout << item->name << " "; });
        cout << endl;
    }

    void printMatrix(Eigen::MatrixXi& matrix) {
        cout << matrix << endl;
    }
#endif

}
