#include <iostream>
#include <string>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <Utils/Alignment.hpp>
#include <ValueObjects/TmdetStruct.hpp>

using namespace std;

namespace Tmdet::Utils::Alignment {

    const int ALIGNMENT_STRIPE = 200; // formerly "alisav"

    vector<string> alignSequences(vector<string> &query, vector<string> &target) {

        vector<string> result;
        Eigen::MatrixXd scores(query.size(), target.size());

        return result;
    }

    vector<string> getChainSequence(
        const Tmdet::ValueObjects::TmdetStruct& tmdetVO, const gemmi::Chain& chain) {

        std::vector<string> sequence;
        auto entityId = chain.residues[0].entity_id;
        for (const auto& entity : tmdetVO.gemmi.entities) {
            if (entity.entity_type == gemmi::EntityType::Polymer
                && entity.name == entityId) {

                sequence = entity.full_sequence;
                break;
            }
        }
        return sequence;
    }

}
