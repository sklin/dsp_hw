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
    //for(int t=0 ; t<sequence.length ; ++t) {
    //    for(int i=0 ; i<hmm->state_num ; ++i) {
    //        cout << alpha[t][i] << "    ";
    //    }
    //    cout << endl;
    //}
    //cout << "===================" << endl;
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
    //for(int t=0 ; t<sequence.length ; ++t) {
    //    for(int i=0 ; i<hmm->state_num ; ++i) {
    //        cout << beta[t][i] << "    ";
    //    }
    //    cout << endl;
    //}
    //cout << "===================" << endl;
}

void gamma_cal(
        const HMM *hmm,
        const Sequence &sequence,
        const double alpha[MAX_SEQ][MAX_STATE],
        const double beta[MAX_SEQ][MAX_STATE],
        double gamma[MAX_SEQ][MAX_STATE],
        double sumGamma[MAX_SEQ][MAX_STATE],
        double sumGammaObserv[MAX_OBSERV][MAX_STATE]
        )
{
    for(int t=0 ; t<sequence.length ; ++t) {
        double denominator = 0.0f;
        for(int i=0 ; i<hmm->state_num ; ++i) {
            denominator += alpha[t][i] * beta[t][i];
        }

        for(int i=0 ; i<hmm->state_num ; ++i) {
            double numerator = alpha[t][i] * beta[t][i];
            if(denominator != 0) {
                gamma[t][i] = numerator / denominator;
                sumGamma[t][i] += gamma[t][i];
                sumGammaObserv[sequence.observ[t]][i] += gamma[t][i];
            }
        }
    }

    //for(int t=0 ; t<sequenceLength ; ++t) {
    //    for(int i=0 ; i<hmm->state_num ; ++i) {
    //        cout << gamma[t][i] << "    ";
    //    }
    //    cout << endl;
    //}
    //cout << "==========================" << endl;
}

void epsilon_cal(
        const HMM *hmm,
        const Sequence &sequence,
        const double alpha[MAX_SEQ][MAX_STATE],
        const double beta[MAX_SEQ][MAX_STATE],
        double epsilon[MAX_SEQ-1][MAX_STATE][MAX_STATE],
        double sumEpsilon[MAX_SEQ][MAX_STATE][MAX_STATE]
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
                if(denominator != 0) {
                    epsilon[t][i][j] = numerator / denominator;
                    sumEpsilon[t][i][j] += epsilon[t][i][j];
                }
            }
        }
        //double prob_sum = 0.0f;
        //for(int i=0 ; i<hmm->state_num ; ++i) {
        //    for(int j=0 ; j<hmm->state_num ; ++j) {
        //        epsilon[t][i][j] = alpha[t][i] * hmm->transition[i][j] * hmm->observation[sequence.observ[t+1]][j] * beta[t+1][j];
        //        prob_sum += epsilon[t][i][j];
        //    }
        //}
        //if(prob_sum!=0) {
        //    for(int i=0 ; i<hmm->state_num ; ++i) {
        //        for(int j=0 ; j<hmm->state_num ; ++j) {
        //            epsilon[t][i][j] /= prob_sum;
        //        }
        //    }
        //}
    }
    //for(int t=0 ; t<sequence.length-1 ; ++t) {
    //    for(int i=0 ; i<hmm->state_num ; ++i) {
    //        for(int j=0 ; j<hmm->state_num ; ++j) {
    //            cout << epsilon[t][i][j] << "    ";
    //        }
    //        cout << endl;
    //    }
    //    cout << "------------------" << endl;
    //}
    //cout << "==================" << endl;
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
        double sumGamma[MAX_SEQ][MAX_STATE] = { {0.0f} };
        double sumGammaObserv[MAX_OBSERV][MAX_STATE] = { {0.0f} };
        double sumEpsilon[MAX_SEQ][MAX_STATE][MAX_STATE] = { { {0.0f} } };

        double updateInitial[MAX_STATE] = {0.0f};
        double updateTransitionNumerator[MAX_STATE][MAX_STATE] = { {0.0f} };
        double updateTransitionDenominator[MAX_STATE] = {0.0f};
        double updateObservationNumerator[MAX_OBSERV][MAX_STATE] = { {0.0f} };
        double updateObservationDenominator[MAX_STATE] = {0.0f};
        for(auto &sequence: sequences) {
            double alpha[MAX_SEQ][MAX_STATE] = { {0.0f} };
            double beta[MAX_SEQ][MAX_STATE] = { {0.0f} };

            double gamma[MAX_SEQ][MAX_STATE] = { {0.0f} };
            double epsilon[MAX_SEQ-1][MAX_STATE][MAX_STATE] = { { {0.0f} } };

            forward_procedure(&hmm, alpha, sequence);
            backward_procedure(&hmm, beta, sequence);
            gamma_cal(&hmm, sequence, alpha, beta, gamma, sumGamma, sumGammaObserv);
            epsilon_cal(&hmm, sequence, alpha, beta, epsilon, sumEpsilon);
        }

        // Update To HMM
        double reverseSeqNum = 1.0f/sequences.size();
        for(int i=0 ; i<hmm.state_num ; ++i)
            hmm.initial[i] = sumGamma[0][i] * reverseSeqNum;

        for(int i=0 ; i<hmm.state_num ; ++i) {
            double denominator = 0.0f;
            for(int t=0 ; t<sequences[0].length-1 ; ++t) { // TODO:
                denominator += sumGamma[t][i];
            }
            for(int j=0 ; j<hmm.state_num ; ++j) {
                double numerator = 0.0f;
                for(int t=0 ; t<sequences[0].length-1 ; ++t) { // TODO:
                    numerator += sumEpsilon[t][i][j];
                }
                if(denominator!=0) {
                    hmm.transition[i][j] = numerator / denominator;
                }
                else {
                    cout << "[Error] Update transition fail." << endl;
                    exit(0);
                }
            }
        }

        for(int i=0 ; i<hmm.state_num ; ++i) {
            double denominator = 0.0f;
            double numerator = 0.0f;
            for(int t=0 ; t<sequences[0].length ; ++t) { // TODO:
                denominator += sumGamma[t][i];
            }
            for(int o=0 ; o<hmm.observ_num ; ++o) {
                numerator = sumGammaObserv[o][i];
                if(denominator!=0) {
                    hmm.observation[o][i] = numerator / denominator;
                }
                else {
                    cout << "[Error] Update observation fail." << endl;
                    exit(0);
                }
            }

        }
    }
    //dumpHMM(stderr, &hmm);
    dumpHMM2File(output_model_name, &hmm);

    return 0;
}
