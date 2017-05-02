#include "hmm.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

#define MAX_MODEL 20

struct Sequence {
    Sequence(string rawData): rawData(rawData), length(rawData.length()) {
        if(rawData.length() > MAX_SEQ) {
            cout << "[Error] Sequence length should not exceed "
                << MAX_SEQ << "!" << endl;
            length = MAX_SEQ;
        }
        for(size_t i=0 ; i<rawData.length() ; ++i)
            observ[i] = rawData[i] - 'A';
        for(size_t i=rawData.length() ; i<MAX_SEQ ; ++i)
            observ[i] = -1;
    }
    string rawData;
    int observ[MAX_SEQ];
    int length;
};


vector<Sequence> load_sequences(const char* filename)
{
    ifstream file(filename);
    if(!file.is_open()) {
        cout << "Sequence file not found!" << endl;
        exit(0);
    }

    vector<Sequence> sequences;
    string line;
    while(getline(file, line)) {
        sequences.emplace_back(line);
    }

    return sequences;
}

double prob(const HMM *hmm, const Sequence &sequence)
{
    double delta[MAX_SEQ][MAX_STATE] = { {0.0f} };
    //int psi[MAX_SEQ][MAX_STATE] = { {0.0f} }; // No need to backtrack the path?

    // initialization
    for(int i=0 ; i<hmm->state_num ; ++i)
        delta[0][i] = hmm->initial[i] * hmm->observation[sequence.observ[0]][i];

    // recursion
    for(int t=0 ; t<sequence.length-1 ; ++t) {
        for(int j=0 ; j<hmm->state_num ; ++j) {
            double maxDeltaAij = 0.0f;
            for(int i=0 ; i<hmm->state_num ; ++i) {
                double deltaAij = delta[t][i] * hmm->transition[i][j];
                if(maxDeltaAij < deltaAij)
                    maxDeltaAij = deltaAij;
            }
            delta[t+1][j] = maxDeltaAij * hmm->observation[sequence.observ[t+1]][j];
        }
    }

    int T = sequence.length-1;
    double maxDelta = 0.0f;
    for(int i=0 ; i<hmm->state_num ; ++i)
        if(maxDelta < delta[T][i])
            maxDelta = delta[T][i];

    return maxDelta;
}


int main(int argc, char **argv)
{
    /**
     * ./test modellist.txt test_data.txt result.txt
     **/
    if(argc < 4) {
        cout << "[Usage] " << argv[0] << " <model_list> <test_data> "
            << "<output>" << endl;
        exit(-1);
    }
    char *model_list = argv[1];
    char *test_data = argv[2];
    char *output = argv[3];

    ofstream file(output);

    HMM models[MAX_MODEL];
    int num = load_models(model_list, models, MAX_MODEL);

    vector<Sequence> sequences(load_sequences(test_data));
    for(auto &sequence: sequences) {
        double maxProb = 0.0f;
        int maxProbIndex = -1;
        for(int n=0 ; n<num ; ++n) {
            double p = prob(&models[n], sequence);
            if(maxProb < p) {
                maxProb = p;
                maxProbIndex = n;
            }
        }
        file << models[maxProbIndex].model_name << " " << maxProb << endl;
        //cout << models[maxProbIndex].model_name << endl;
        //fprintf(fp, "%s %f\n", models[maxProbIndex].model_name, maxProb);
    }

    return 0;
}
