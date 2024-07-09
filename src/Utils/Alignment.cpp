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

    void alignSequences(gemmi::Chain &query, vector<string> &target) {

        int alignmentStripe = 200; // formerly "alisav"
        int maxgap = 300;
        // number of residues based on atom lines:
        int queryLength = query.residues.size();
        // number of residues based on SEQRES/entity_poly info:
        int targetLength = target.size();

        Eigen::MatrixXd scores(targetLength, queryLength);
        Eigen::MatrixXd iPath(targetLength, queryLength);
        Eigen::MatrixXd jPath(targetLength, queryLength);

        if (abs(queryLength - targetLength) > 200) {
            alignmentStripe = abs(queryLength - targetLength) + 200;
        }
        if (abs(queryLength - targetLength) > 100) {
            maxgap = 1000;
        }
        vector<int> gaps;
        gaps.resize(queryLength);
        // gaps init
        initGaps(query, gaps);

        // Calculating scores, iPath and jPath
        const int GAPOPEN = 2;
        int maxScore, i, j;
        for (i = 0; i < targetLength; i++) {
            for (j = (i-alignmentStripe < 0 ? 0 : i-alignmentStripe);
                (j < queryLength && j < i+alignmentStripe); j++) {

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
                if (j>0) {
                    for (int k = i-2; (k>=0 && k> i-maxgap); k--) {
                        int gap = gaps[j-1] * (i-k) + GAPOPEN;
                        if (scores(k, j-1) - gap >= maxScore) {
                            maxScore = scores(k, j-1) - gap;
                            iPath(i, j) = k;
                        }
                    }
                }
                scores(i, j) = (target[i] == query.residues[j].name ? 10 : 0) + maxScore;
            }
        }
        maxScore = scores(targetLength-1, queryLength-1);
        list<string> targetList(target.begin(), target.end());
        list<gemmi::Residue> queryList(query.residues.begin(), query.residues.end());

        // trace back paths?
        int nuti = targetLength - 1;
        for (int i = nuti; i >= 0; i--) {
            if (maxScore <= scores(i, queryLength - 1)) {
                maxScore = scores(i, queryLength - 1);
                nuti = i;
            }
        }
        int nutj = queryLength - 1;
        for (int j = nutj; j >= 0; j--) {
            if (maxScore <= scores(targetLength - 1, j)) {
                maxScore = scores(targetLength - 1, j);
                nutj = j;
            }
        }

        // update target sequence
        for (int j = nutj+1; j < queryLength; j++) {
            targetList.emplace_back(query.residues[j].name);
        }

        // inserts
        int UTI = nuti;
        int UTJ = nutj;

        while (UTI >= 0 && UTJ >= 0) {
            auto currentResidue = query.residues[UTJ];
            if (target[UTI] != currentResidue.name) {
                cerr << "Conflict in residue name: chain "
                    << query.name << " at position " << currentResidue.label_seq.value
                    << target[UTI] << " vs " << currentResidue.name << endl;
                target[UTI] = currentResidue.name;
            }
            nuti = iPath(UTI, UTJ);
            nutj = jPath(UTI, UTJ);
            for (int j = UTJ-1; j > nutj; j--) {
                // insert before
                // TODO !!!!!!!!!!!!!!!!!!!!!
            }
        }


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
        // std::sort(chain.residues.begin(), chain.residues.end(), compareResidues);
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
            int residueDistance = current->label_seq.value - previous->label_seq.value;
            if (residueDistance > 1 && residueDistance < 500) {
                atomDistance = 20;
            }
            if (atomDistance > 18.3) {
                if (index <= GAP) {
                    gaps[index-1] = 2;
                } else {
                    gaps[index-1] = 0;
                }
                if (index >= gaps.size() - GAP) {
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

}
