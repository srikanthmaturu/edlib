#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <ctime>
#include <string>
#include <climits>
#include <queue>
#include <fstream>
#include <iostream>
#include <edlib.h>
#include <cmath>
#include "edlib.h"

using namespace std;

int readFastaSequences(const char* path, vector< vector<char> >* seqs);

void printAlignment(const char* query, const char* target,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const EdlibAlignMode modeCode);

// For debugging
void printSeq(const vector<char> &seq) {
    for (int i = 0; i < (int) seq.size(); i++)
        printf("%d ", seq[i]);
    printf("\n");
}

EdlibEqualityPair additionalEqualities[4] = {{'B','N'},{'Z','Q'}, {'x','A'}, {'X','A'}};

uint64_t fastPercentIdentity(string s1, string s2, uint64_t percentIdentityThreshold) {
    uint64_t s1_size = s1.size(), s2_size = s2.size();
    double den = (s1_size > s2_size) ? s1_size : s2_size;
    double num = ((s1_size > s2_size)) ? s2_size : s1_size;
    //double pi = (percentIdentityThreshold * 1.0)/100.0;

    //uint64_t thEd = ceil((num - pi*den) / (double)pi);

    EdlibAlignConfig edlibConfig = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, additionalEqualities, 4);
    EdlibAlignResult ed_result = edlibAlign(s1.c_str(), s1.size(), s2.c_str(), s2.size(), edlibConfig);

    uint64_t matches = 0;
    for(int64_t i = 0; i < ed_result.alignmentLength; i++) {
        if(ed_result.alignment[i] == EDLIB_EDOP_MATCH) {
            matches++;
        }
    }
    auto p =  round((matches * 100.0)/(double)ed_result.alignmentLength);
    edlibFreeAlignResult(ed_result);
    return p;
}

uint64_t fastInfixPercentIdentity(string s1, string s2) {
    uint64_t s1_size = s1.size(), s2_size = s2.size();
    double den = (s1_size > s2_size) ? s1_size : s2_size;
    double num = ((s1_size > s2_size)) ? s2_size : s1_size;

    EdlibAlignConfig edlibConfig = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 4);
    EdlibAlignResult ed_result = edlibAlign(s1.c_str(), s1.size(), s2.c_str(), s2.size(), edlibConfig);

    uint64_t matches = 0;
    for(int64_t i = 0; i < ed_result.alignmentLength; i++) {
        if(ed_result.alignment[i] == EDLIB_EDOP_MATCH) {
            matches++;
        }
    }
    auto p =  round((matches * 100.0)/(double)ed_result.alignmentLength);
    edlibFreeAlignResult(ed_result);
    return p;
}

string alignSequences(vector<char> firstSequence, vector<char> secondSequence, EdlibAlignConfig edlibAlignConfig){
    //cout << string{firstSequence.begin(), firstSequence.end()} << " \n" << string{secondSequence.begin(), secondSequence.end()} << endl;
    EdlibAlignResult alignResult = edlibAlign(firstSequence.data(), firstSequence.size(), secondSequence.data(), secondSequence.size(), edlibAlignConfig);
    int tIdx = alignResult.endLocations[0], qIdx = -1;
    for (int i = 0; i < alignResult.alignmentLength; i++) {
        if (alignResult.alignment[i] != EDLIB_EDOP_INSERT)
            tIdx--;
    }
    //cout << alignResult.startLocations[0] << "\t" << tIdx << "\t" <<endl;
    string alignSequence;
    if(tIdx >= 0){
        alignSequence = {secondSequence.begin(), secondSequence.begin() + tIdx + 1};
    }
    auto alignment = alignResult.alignment;
    for (int i = 0; i < alignResult.alignmentLength; i++) {
        switch(alignment[i]){
            case EDLIB_EDOP_DELETE:
                alignSequence = alignSequence + secondSequence[++tIdx];
                //cout << "DELETE" << "\t";
                break;
            case EDLIB_EDOP_INSERT:
                alignSequence = alignSequence + firstSequence[++qIdx];
                //cout << "INSERT" << "\t";
                break;
            case EDLIB_EDOP_MATCH:
                alignSequence = alignSequence + secondSequence[++tIdx];
                ++qIdx;
                //cout << "MATCH" << "\t";
                break;
            case EDLIB_EDOP_MISMATCH:
                alignSequence = alignSequence + secondSequence[++tIdx];
                ++qIdx;
                //cout << "MISMATCH" << "\t";
        }
    }
    cout << endl;
    if(tIdx < secondSequence.size()){
        alignSequence = alignSequence + string{secondSequence.begin() + tIdx + 1, secondSequence.end()};
    }
    printAlignment(firstSequence.data(), secondSequence.data(),  alignResult.alignment, alignResult.alignmentLength,  *(alignResult.endLocations), EDLIB_MODE_HW);
    cout << "Edit distance: " << alignResult.editDistance << endl;
    edlibFreeAlignResult(alignResult);
    cout << "Percent identity: " << fastPercentIdentity(firstSequence.data(), secondSequence.data(), 0) << endl;
    cout << "Infix Percent identity: " << fastInfixPercentIdentity(firstSequence.data(), secondSequence.data()) << endl;
    return alignSequence;
}

int main(int argc, char* const argv[]){
    if(argc < 3){
        cout << "Usage ./aligner2 [sequences_file] [output_file]" << endl;
        return 1;
    }
    EdlibAlignConfig edlibConfig = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0);

    int readResult;
    // Read queries
    char* sequenceFilepath = argv[1];
    vector< vector<char> >* sequences = new vector< vector<char> >();
    printf("Reading sequences... \n");
    readResult = readFastaSequences(sequenceFilepath, sequences);
    cout << "Number of sequences: " << (sequences)->size() << endl;
    if (readResult) {
        printf("Error: There is no file with name %s\n", sequenceFilepath);
        delete sequences;
        return 1;
    }
    ofstream results(argv[2], ofstream::out);

    for(uint64_t i = 0; i < sequences->size(); i++){
        for(uint64_t j = i + 1; j < sequences->size(); j++){
            results << ">Alignment of sequences " << i << " and " << j << endl;
            cout << ">Alignment of sequences " << i << " and " << j << endl;
            string alignedSequence = alignSequences((*sequences)[i], (*sequences)[j], edlibConfig);
            results << alignedSequence << endl;
        }
    }
}

/** Reads sequences from fasta file.
 * @param [in] path Path to fasta file containing sequences.
 * @param [out] seqs Sequences will be stored here, each sequence as vector of letters.
 * @return 0 if all ok, positive number otherwise.
 */
int readFastaSequences(const char* path, vector< vector<char> >* seqs) {
    seqs->clear();

    FILE* file = fopen(path, "r");
    if (file == 0)
        return 1;

    bool inHeader = false;
    bool inSequence = false;
    const int buffSize = 4096;
    char buffer[buffSize];
    while (!feof(file)) {
        int read = fread(buffer, sizeof(char), buffSize, file);
        for (int i = 0; i < read; ++i) {
            char c = buffer[i];
            if (inHeader) { // I do nothing if in header
                if (c == '\n')
                    inHeader = false;
            } else {
                if (c == '>') {
                    inHeader = true;
                    inSequence = false;
                } else {
                    if (c == '\r' || c == '\n')
                        continue;
                    // If starting new sequence, initialize it.
                    if (inSequence == false) {
                        inSequence = true;
                        seqs->push_back(vector<char>());
                    }
                    seqs->back().push_back(c);
                }
            }
        }
    }

    fclose(file);
    return 0;
}


void printAlignment(const char* query, const char* target,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const EdlibAlignMode modeCode) {
    int tIdx = -1;
    int qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = position;
        for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] != EDLIB_EDOP_INSERT)
                tIdx--;
        }
    }
    for (int start = 0; start < alignmentLength; start += 50) {
        // target
        printf("T: ");
        int startTIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_INSERT)
                printf("-");
            else
                printf("%c", target[++tIdx]);
            if (j == start)
                startTIdx = tIdx;
        }
        printf(" (%d - %d)\n", max(startTIdx, 0), tIdx);

        // match / mismatch
        printf("   ");
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
        }
        printf("\n");

        // query
        printf("Q: ");
        int startQIdx = qIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_DELETE)
                printf("-");
            else
                printf("%c", query[++qIdx]);
            if (j == start)
                startQIdx = qIdx;
        }
        printf(" (%d - %d)\n\n", max(startQIdx, 0), qIdx);
    }
}
