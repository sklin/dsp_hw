#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <string>
#include <map>
#include <algorithm> // for max_element

#include <stdio.h>

#include "Ngram.h"

#define DEBUG_MESSAGE

using namespace std;

void printUsage()
{
    cout << "[Usage] ./disambig -text $file -map $map -lm $LM -order $order" << endl;
}

#define Big5 unsigned

Big5 CharToBig5(char b0, char b1)
{
    Big5 _b0 = (unsigned) b0 & 0xff;
    Big5 _b1 = (unsigned) b1 & 0xff;
    Big5 result = (_b0 << 8) | _b1;
    //printf("0x%hX,\t0x%hX=\t%hX\n", _b0 << 8, _b1, result);
    return result;
}

void Big5ToChar(Big5 w, char *str)
{
    str[0] = w >> 8;
    str[1] = w & 0xff;
    str[2] = '\0';
}

/*
void Big5ToChar(Big5 w, char *b0, char *b1)
{
    *b0 = (w >> 8) & 0xff;
    *b1 = w & 0xff;
}
*/

Big5 StringToBig5(const string& str)
{
    if(str.length() != 2) {
        cout << "Detect non-Big5 word!" << endl;
        exit(-1);
    }
    return CharToBig5(str[0], str[1]);
}

class MyDisambig {
public:
	MyDisambig(const char*, const char*, const char*, int);
	~MyDisambig();
private:
	VocabIndex getIndex(const char *word) { return voc.getIndex(word); }
	VocabIndex getIndex(char b0, char b1) {
		//char tmp[3] = { b0, b1, '\0' };
		char tmp[3];
        tmp[0] = b0; tmp[1] = b1; tmp[2] = '\0';
		return voc.getIndex(tmp);
	}
    VocabIndex getIndex(Big5 w) {
        char tmp[3];
        //Big5ToChar(w, &tmp[0], &tmp[1]);
        tmp[0] = w >> 8;
        tmp[1] = w & 0xff;
        tmp[2] = '\0';
        return voc.getIndex(tmp);
    }

    double getUnigramProb(const char *w1);
    double getBigramProb(const char *w1, const char *w2);
    double getTrigramProb(const char *w1, const char *w2, const char *w3);
    vector<Big5> Viterbi(vector<Big5>& seq);

	Ngram *lm;
    Vocab voc;
    map<Big5, vector<Big5> > mapping;
    vector<string> content;
	int order;
};

MyDisambig::MyDisambig(const char *textFilename, const char *mapFilename, const char *lmFilename, int order)
	: order(order)
{
    // Load text file
    fstream inFile(textFilename);
    if(!inFile.is_open()) {
        cout << "Cannot open " << textFilename << endl;
        exit(-1);
    }
    while(!inFile.eof()) {
        string line;
        getline(inFile, line);
        if(line.length() > 0)
            content.push_back(line);
    }
    inFile.close();
#ifdef DEBUG_MESSAGE
    cout << "content.size(): " << content.size() << endl;
#endif

    // Load Mapping
    fstream mapFile(mapFilename);
    if(!mapFile.is_open()) {
        cout << "Cannot open " << mapFilename << endl;
        exit(-1);
    }
    int cnt = 0;
    while(!mapFile.eof()) {
        string line;
        getline(mapFile, line);
        if(line.length() == 0)
            break;
        stringstream ss(line);

        string key, value;
        vector<Big5> values;
        ss >> key;
        while(ss >> value) {
            values.push_back(StringToBig5(value));
        }
        mapping[StringToBig5(key)] = values;
        cnt++;
    }
#ifdef DEBUG_MESSAGE
    cout << "mapping.size(): " << mapping.size() << endl;
    cout << "cnt: " << cnt << endl;
#endif

    // Load LM
    lm = new Ngram(voc, order);
    File lmFile(lmFilename, "r");
    lm->read(lmFile);
    lmFile.close();
}

MyDisambig::~MyDisambig()
{
    delete lm;
}

vector<Big5> MyDisambig::Viterbi(vector<Big5> &seq)
{
    vector<Big5> out;
    int seqSize = seq.size();

    // Allocate
    // Because #states of each observ. is not the same,
    //  I use std::vector to store delta values.
    vector<vector<double> > delta;
    for(int i=0 ; i<seqSize ; ++i) {
        delta[i].resize(mapping[seq[i]].size(), 0.0f);
    }

    // Initialize, t=0
    Big5 seq0 = seq[0];
    VocabIndex wid = getIndex(seq0);
    for(int i=0 ; i<mapping[seq0].size() ; ++i) { // mapping[seq0].size() == delta[0].size()
        char w[3];
        Big5ToChar(mapping[seq0][i], w);
        delta[0][i] = getUnigramProb(w);
    }

    // Recursion, t=1
    for(int t=1 ; t<seqSize ; ++t) {
        for(int i=0 ; i<mapping[seq[t]].size() ; ++i) { // mapping[seq[t]].size() == delta[t].size()
            // find max getBigramProb(seq[) j for j=0~delta[t-1].size()-1
            double maxProb = 0.0f;
            for(int j=0 ; j<mapping[seq[t-1]].size() ; ++j) {
                char w1[3], w2[3];
                Big5ToChar(mapping[seq[t]][i], w1);
                Big5ToChar(mapping[seq[t-1]][j], w2);
                double prob = getBigramProb(w1, w2);
                if(prob > maxProb)
                    maxProb = prob;
            }
            delta[t][i] = maxProb;
        }
    }

    double finalProb = *max_element(delta[seqSize-1].begin(), delta[seqSize-1].end());
}

// Get P(W1) -- unigram
double MyDisambig::getUnigramProb(const char *w1)
{
    VocabIndex wid1 = voc.getIndex(w1);

    if(wid1 == Vocab_None)  //OOV
        wid1 = voc.getIndex(Vocab_Unknown);

    VocabIndex context[] = { Vocab_None };
    return lm->wordProb( wid1, context);
}

// Get P(W2 | W1) -- bigram
double MyDisambig::getBigramProb(const char *w1, const char *w2)
{
    VocabIndex wid1 = voc.getIndex(w1);
    VocabIndex wid2 = voc.getIndex(w2);

    if(wid1 == Vocab_None)  //OOV
        wid1 = voc.getIndex(Vocab_Unknown);
    if(wid2 == Vocab_None)  //OOV
        wid2 = voc.getIndex(Vocab_Unknown);

    VocabIndex context[] = { wid1, Vocab_None };
    return lm->wordProb( wid2, context);
}

// Get P(W3 | W1, W2) -- trigram
double MyDisambig::getTrigramProb(const char *w1, const char *w2, const char *w3) 
{
    VocabIndex wid1 = voc.getIndex(w1);
    VocabIndex wid2 = voc.getIndex(w2);
    VocabIndex wid3 = voc.getIndex(w3);

    if(wid1 == Vocab_None)  //OOV
        wid1 = voc.getIndex(Vocab_Unknown);
    if(wid2 == Vocab_None)  //OOV
        wid2 = voc.getIndex(Vocab_Unknown);
    if(wid3 == Vocab_None)  //OOV
        wid3 = voc.getIndex(Vocab_Unknown);

    VocabIndex context[] = { wid2, wid1, Vocab_None };
    return lm->wordProb( wid3, context );
}

int main(int argc, char **argv)
{
    //          0     1     2    3    4   5   6      7      8
    // ./disambig -text $file -map $map -lm $LM -order $order
    char *textFilename = NULL, *mapFilename = NULL, *lmFilename = NULL;
    int order = 2;

    int argCnt = 1;
    while(argCnt+1 < argc) {
        if(strcmp(argv[argCnt], "-text") == 0) {
            textFilename = argv[argCnt+1];
        }
        else if(strcmp(argv[argCnt], "-map") == 0) {
            mapFilename = argv[argCnt+1];
        }
        else if(strcmp(argv[argCnt], "-lm") == 0) {
            lmFilename = argv[argCnt+1];
        }
        else if(strcmp(argv[argCnt], "-order") == 0) {
            order = atoi(argv[argCnt+1]);
        }
        else {
            printUsage();
            exit(-1);
        }
        argCnt += 2;
    }

    if(!textFilename || !mapFilename || !lmFilename) {
        printUsage();
        exit(-1);
    }

#ifdef DEBUG_MESSAGE
    cout << "file: " << textFilename << endl;
    cout << "map: " << mapFilename << endl;
    cout << "lm: " << lmFilename << endl;
    cout << "order: " << order << endl;
#endif
    
    MyDisambig myDisamdig(textFilename, mapFilename, lmFilename, order);
    return 0;
}
