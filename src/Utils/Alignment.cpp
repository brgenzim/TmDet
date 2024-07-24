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
    static gemmi::Residue* copyResidue(const gemmi::Residue& other);
    static void initGaps(gemmi::Chain &chain, vector<int> &gaps);

#ifdef __ALIGNMENT_DBG
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

    string currentCode;

    static void printMatrix(Eigen::MatrixXi& matrix);
    static void writeMatrix(FILE* file, Eigen::MatrixXi *mx, int rows, int cols);
    static void writeDebugState(alignState* state);
    void printResidues(const list<gemmi::Residue> &residues);
    void printResidues(const vector<gemmi::Residue> &residues);
    void printResidues(const vector<gemmi::Residue*> &residues);
    void printResidues(const list<gemmi::Residue*> &residues);

    #define MATRIX_DIR "/tmp/pdb_map"
#endif

    void alignSequences(gemmi::Chain &chain, vector<string> sequence) {

#ifdef __ALIGNMENT_DBG
        cout << "WARNING: only chain G will be aligned in DEBUG mode" << endl;
        if (chain.name != "G") {
            return;
        }
#endif

        vector<gemmi::Residue *> seqResidues;
        int seqNum = 1;
        for (auto &residueName : sequence) {
            auto residue = createResidue(-9999, seqNum, residueName, chain.name);
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
        // TODO: remove later
        // This entry code is hardwired since the function does not have any extra parameter for it.
        alignState state;
        state.code = currentCode;
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

        maxScore = scores(seqResiduesLength-1, chainLength-1);
        // Create lists from sequence and atom lines
        list<gemmi::Residue *> seqResiduesList(seqResidues.begin(), seqResidues.end());
        // residues from atom lines
        list<gemmi::Residue> chainList;
        for (auto it = chain.residues.begin(); it < chain.residues.end(); it++) {
            auto copy = it->empty_copy();
            copy.atoms = it->atoms;
            chainList.push_back(copy);
        }

        // trace back paths?
        int UTI = seqResiduesLength;
        int UTJ = chainLength;
        int nuti = seqResiduesLength - 1;
        for (int i = nuti; i >= 0; i--) {
#ifdef __ALIGNMENT_DBG
            cout << "scores(" << i << ", " << chainLength-1 << "): " << scores(i, chainLength-1);
#endif
            if (maxScore <= scores(i, chainLength - 1)) {
                maxScore = scores(i, chainLength - 1);
                nuti = i;
#ifdef __ALIGNMENT_DBG
                cout << " nuti: " << i << " maxScore: " << maxScore;
#endif
            }
#ifdef __ALIGNMENT_DBG
            cout << endl;
#endif
        }
        int nutj = chainLength - 1;
        for (int j = nutj; j >= 0; j--) {
            if (maxScore <= scores(seqResiduesLength - 1, j)) {
                maxScore = scores(seqResiduesLength - 1, j);
                nutj = j;
                nuti = seqResiduesLength - 1;
            }
        }

        // update seqResidues sequence
        for (int j = nutj+1; j < chainLength; j++) {
            auto residueFromAtomLine = copyResidue(chain.residues[j]);
            seqResiduesList.push_back(residueFromAtomLine);
        }

        // inserts
        UTI = nuti;
        UTJ = nutj;
#ifdef __ALIGNMENT_DBG
        std::cout << "UTI: " << UTI << " UTJ: " << UTJ << " nuti: " << nuti << " nutj: " << nutj << endl;
#endif

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
                auto residueFromAtomLine = copyResidue(chain.residues[j]);
                // lambda expression for find_if call
                auto equals = [&seqResidues, &UTI](gemmi::Residue* item) {
                    return item->seqid.num == seqResidues[UTI]->seqid.num;
                };
                auto pos = find_if(seqResiduesList.begin(), seqResiduesList.end(), equals);
                // insert pointer before position
                seqResiduesList.insert(pos, residueFromAtomLine);
            }
            UTI = nuti;
            UTJ = nutj;
        }

#ifdef __ALIGNMENT_DBG
        printf("Chain %s: UTI: %d, UTJ: %d, nuti: %d, nutj: %d\n",
            chain.name.c_str(), UTI, UTJ, nuti, nutj);

        cout << "Chain " << chain.name << ": seqResiduesList: ";
        printResidues(seqResiduesList);
        cout << "Chain " << chain.name << ": atom residues: ";
        printResidues(chainList);
#endif

        // Filling gaps in residue list from atom lines

        // auto seqResidue = seqResiduesList.begin();
        // auto chainResidue = chainList.begin();
        // while (seqResidue != seqResiduesList.end()) {

        //     if (chainResidue == chainList.end()) {
        //         chainList.emplace_back(*seqResidue);
        //         ++seqResidue;
        //         continue;
        //     }

        //     auto seqLabelValue = seqResidue->label_seq.value;
        //     auto chainLabelValue = chainResidue->label_seq.value;

        //     if (seqLabelValue < chainLabelValue) {
        //         // insert does not invalidate iterators
        //         chainList.insert(chainResidue, *seqResidue);
        //         ++seqResidue;
        //     } else if (seqLabelValue == chainLabelValue) {
        //         ++seqResidue;
        //         ++chainResidue;
        //     } else { // chain label < seq label
        //         ++chainResidue;
        //         chainList.insert(chainResidue, *seqResidue);
        //         --chainResidue;
        //     }
        // }

        // Update residues of chain:
        // rescue data from chainList pointers - since most of its items
        // still reference items in chain.residues and that will be cleared
        chain.residues.clear();
        for (auto residue : chainList) {
            auto copy = residue.empty_copy();
            copy.atoms = residue.atoms;
            chain.residues.emplace_back(copy);
        }
        // clean ups
        auto free = [](gemmi::Residue* item) { delete item; };
        for_each(seqResiduesList.begin(), seqResiduesList.end(), free);
        seqResiduesList.resize(0);
        seqResidues.resize(0);
    }

    gemmi::Residue* createResidue(int seqNum, int labelSeqNum, string name, string chainName) {
        auto residue = new gemmi::Residue();
        residue->seqid.num.value = seqNum;
        residue->label_seq.value = labelSeqNum;
        residue->name = name;
        residue->subchain = chainName;
        return residue;
    }

    gemmi::Residue* copyResidue(const gemmi::Residue& other) {
        auto residue = new gemmi::Residue();
        residue->seqid = other.seqid;
        residue->label_seq = other.label_seq;
        residue->name = other.name;
        residue->subchain = other.subchain;

        residue->entity_id = other.entity_id;
        residue->entity_type = other.entity_type;
        residue->atoms = other.atoms;
        residue->flag = other.flag;
        residue->group_idx = other.group_idx;
        residue->het_flag = other.het_flag;
        residue->segment = other.segment;
        residue->sifts_unp = other.sifts_unp;

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

#ifdef __ALIGNMENT_DBG
            currentCode = tmdetVO.gemmi.name;
            for (unsigned int i = 0; i < currentCode.size(); ++i) {
                currentCode[i] = tolower(currentCode[i]);
            }
#endif

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

#ifdef __ALIGNMENT_DBG
    void _printResidue(gemmi::Residue item) {
        cout << item.name << "[" << item.seqid.num.value << "] ";
    }

    void _printResidue(gemmi::Residue* item) {
        cout << item->name << "[" << item->seqid.num.value << "] ";
    }

    void printResidues(const vector<gemmi::Residue> &residues) {
        auto start = residues.begin();
        auto end = residues.end();
        while (start != end) {
            _printResidue(*start);
            ++start;
        }
        cout << endl;
    }

    void printResidues(const list<gemmi::Residue> &residues) {
        auto start = residues.begin();
        auto end = residues.end();
        while (start != end) {
            _printResidue(*start);
            ++start;
        }
        cout << endl;
    }

    void printResidues(const vector<gemmi::Residue*> &residues) {
        auto start = residues.begin();
        auto end = residues.end();
        while (start != end) {
            _printResidue(*start);
            ++start;
        }
        cout << endl;
    }

    void printResidues(const list<gemmi::Residue*> &residues) {
        auto start = residues.begin();
        auto end = residues.end();
        while (start != end) {
            _printResidue(*start);
            ++start;
        }
        cout << endl;
    }

    void printMatrix(Eigen::MatrixXi& matrix) {
        cout << matrix << endl;
    }

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

        if (currentCode != "7ec3" || state->chain != 'G') {
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

}
