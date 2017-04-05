#include "hmm.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

/**
 * alpha, beta, gamma, epsilon
 **/


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

void forward_procedure(const HMM* hmm, double alpha[MAX_SEQ][MAX_STATE], const Sequence &sequence)
{
    /**
     * alpha[t][i] denote alpha_t(i),
     *  where `i` is state 0 ~ state hmm->state_num
     *  and `t` is the time for 0 ~ sequence.length-1
     **/

    for(int i=0 ; i<hmm->state_num ; ++i) {
        alpha[0][i] = hmm->initial[i] * hmm->observation[sequence.observ[0]][i];
    }

    for(int t=1 ; t<sequence.length ; ++t) {
        for(int i=0 ; i<hmm->state_num ; ++i) {
            alpha[t][i] = 0.0f;
            for(int j=0 ; j<hmm->state_num ; ++j) {
                alpha[t][i] += alpha[t-1][j] * hmm->transition[j][i];
            }
            alpha[t][i] *= hmm->observation[sequence.observ[t]][i];
        }
    }
}

void backward_procedure(const HMM* hmm, double beta[MAX_SEQ][MAX_STATE], const Sequence &sequence)
{
    /**
     * beta[t][i] denote beta_t(i), where 1 <= i <= N
     **/
    for(int i=0 ; i<hmm->state_num ; ++i) {
        beta[sequence.length-1][i] = 1;
    }

    for(int t=sequence.length-2 ; t>=0 ; --t) {
        for(int i=0 ; i<hmm->state_num ; ++i) {
            beta[t][i] = 0.0f;
            for(int j=0 ; j<hmm->state_num ; ++j) {
                beta[t][i] += hmm->transition[i][j] 
                    * hmm->observation[sequence.observ[t+1]][j] 
                    * beta[t+1][j];
            }
        }
    }
}

void gamma_cal(
        const HMM *hmm,
        int sequenceLength,
        const double alpha[MAX_SEQ][MAX_STATE],
        const double beta[MAX_SEQ][MAX_STATE],
        double gamma[MAX_SEQ][MAX_STATE]
        )
{
    for(int t=0 ; t<sequenceLength ; ++t) {
        double denominator = 0.0f;
        for(int j=0 ; j<hmm->state_num ; ++j) {
            denominator += alpha[t][j] * beta[t][j];
        }

        for(int i=0 ; i<hmm->state_num ; ++i) {
            double numerator = alpha[t][i] * beta[t][i];
            if(denominator == 0) { // TODO:!
                //cout << "gamma_cal: demonicator=0, numerator=" << numerator << endl;
                //exit(-1);
            }
            else {
                gamma[t][i] = numerator / denominator;
            }
        }

    }
}

void epsilon_cal(
        const HMM *hmm,
        const Sequence &sequence,
        const double alpha[MAX_SEQ][MAX_STATE],
        const double beta[MAX_SEQ][MAX_STATE],
        double epsilon[MAX_SEQ-1][MAX_STATE][MAX_STATE]
        )
{
    for(int t=0 ; t<sequence.length-1 ; ++t) {
        double denominator = 0.0f;
        for(int i=0 ; i<hmm->state_num ; ++i) {
            for(int j=0 ; j<hmm->state_num ; ++j) {
                denominator += alpha[t][i] * hmm->transition[i][j] * hmm->observation[sequence.observ[t+1]][j] * beta[t+1][j];
            }
        }

        for(int i=0 ; i<hmm->state_num ; ++i) {
            for(int j=0 ; j<hmm->state_num ; ++j) {
                double numerator = alpha[t][i] * hmm->transition[i][j] * hmm->observation[sequence.observ[t+1]][j] * beta[t+1][j];
                if(denominator == 0) { // TODO:!
                    //cout << "epsilon_cal: demonicator=0, numerator=" << numerator << endl;
                    //exit(-1);
                }
                else {
                    epsilon[t][i][j] = numerator / denominator;
                }
            }
        }
    }
}

void update(
        HMM *hmm,
        const double gamma[MAX_SEQ][MAX_STATE],
        const double epsilon[MAX_SEQ-1][MAX_STATE][MAX_STATE],
        const Sequence &sequence
        )
{
    // hmm->initial
    for(int i=0 ; i<hmm->state_num ; ++i)
        hmm->initial[i] = gamma[0][i];

    // hmm->transition
    for(int i=0 ; i<hmm->state_num ; ++i) {
        double denominator = 0.0f;
        for(int t=0 ; t<sequence.length-1 ; ++t) {
            denominator += gamma[t][i];
        }
        for(int j=0 ; j<hmm->state_num ; ++j) {
            double numerator = 0.0f;
            for(int t=0 ; t<sequence.length-1 ; ++t) {
                numerator += epsilon[t][i][j];
            }
            if(denominator == 0) { // TODO:!
                //cout << "train_hmm: demonicator=0, numerator=" << numerator << endl;
                //exit(-1);
            }
            else {
                hmm->transition[i][j] = numerator / denominator;
            }
        }
    }

    // hmm->observation
    for(int i=0 ; i<hmm->state_num ; ++i) {
        double numerators[MAX_OBSERV] = { 0.0f };
        double denominator = 0.0f;
        for(int t=0 ; t<sequence.length ; ++t) { // TODO: t: sequence.length & T: MAX_SEQ ?
            numerators[sequence.observ[t]] += gamma[t][i];
            denominator += gamma[t][i];
        }
        if(denominator == 0) { // TODO:!
            //cout << "train_hmm: demonicator=0" << endl;
            //exit(-1);
        }
        else {
            for(int k=0 ; k<MAX_OBSERV ; ++k) {
                hmm->observation[k][i] = numerators[k] / denominator;
            }
        }
    }
    //for(int k=0 ; k<hmm->observ_num ; ++k) {
    //    for(int j=0 ; j<hmm->state_num ; ++j) {
    //        double numerators[MAX_STATE] = { 0.0f };
    //        double denominator = 0.0f;
    //        for(int t=0 ; t<sequence.length ; ++t) {
    //            denominator += gamma[t][j];
    //            numerators[sequence.observ[t]] += gamma[t][j];
    //        }
    //        hmm->observation[k][j] = numerator / denominator;
    //    }
    //}
}

int main(int argc, char **argv)
{
    /**
     * ./train iteration model_init.txt seq_model_*.txt model_*.txt
     **/
    if(argc < 5) {
        cout << "[Usage] " << argv[0] << " <iteration> <model_init> "
            << "<input_model_name> <output_model_name>" << endl;
        exit(-1);
    }
    int iteration = atoi(argv[1]);
    char *model_init = argv[2];
    char *input_model_name = argv[3];
    char *output_model_name = argv[4];

	HMM hmm;
	loadHMM(&hmm, model_init);
    vector<Sequence> sequences(load_sequences(input_model_name));

    for(int n=0 ; n<iteration ; ++n) {
        for(auto &sequence: sequences) {
            double alpha[MAX_SEQ][MAX_STATE] = { {0.0f} };
            double beta[MAX_SEQ][MAX_STATE] = { {0.0f} };
            double gamma[MAX_SEQ][MAX_STATE] = { {0.0f} };
            double epsilon[MAX_SEQ-1][MAX_STATE][MAX_STATE] = { { {0.0f} } };
            forward_procedure(&hmm, alpha, sequence);
            backward_procedure(&hmm, beta, sequence);
            gamma_cal(&hmm, sequence.length, alpha, beta, gamma);
            epsilon_cal(&hmm, sequence, alpha, beta, epsilon);

            for(int i=0 ; i<sequence.length ; ++i) {
                for(int j=0 ; j<hmm.state_num ; ++j)
                    cout << gamma[i][j] << "    ";
                cout << endl;
            }

            update(&hmm, gamma, epsilon, sequence);
        }
    }
    dumpHMM(stderr, &hmm);

    return 0;
}
