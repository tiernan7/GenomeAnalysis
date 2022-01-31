#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <chrono>
using namespace std;

//A struct for working with matches
//Helper functions for string manipulation

char getComplement(char base){ //returns the reverse complement of a base
    if (base == 'A'){
        return 'T';
    }
    if (base == 'C'){
       return 'G';
    }
    if (base == 'G'){
       return 'C';
    }
    if (base == 'T'){
       return 'A';
    }
    else{
        return 'N';
    }
}


string * generateReverseComplement(string * strand){ //allocates the reverse complement of a strand to the heap
    string * outstr = new string();
    *outstr = "";
    string str = *strand;
    for (int i = 0; i < (*strand).size(); i++){
        (*outstr).push_back(getComplement(str[i]));
    }
    reverse((*outstr).begin(), (*outstr).end());
    return outstr;
}

int countBases(string * seq, char base){ //counts the number of occurences of a character in a string
    int count = 0;
    for (int i = 0; i < (*seq).size(); i++){
        if ((*seq)[i] == base){
            count++;
        }
    }
    return count;
}

string * readFASTA(string fileName){ //reads a standard fasta file, allocates it on the heap, and returs a pointer to its location
    string firstLine;
    string *returnString = new string;
    ifstream file;
    string buffer;
    file.open(fileName);
    if (file.is_open()){
        getline(file,firstLine);
        while (getline(file,buffer)){
            (*returnString).append(buffer);
        }
    }
    transform((*returnString).begin(), (*returnString).end(), (*returnString).begin(), ::toupper);
    return (returnString);
}
//Forward declarations for global variables
int compare (const void * a, const void * b);
ofstream outfile;


//string * sequence = readFASTA("test1.txt");
//string * sequence2 = readFASTA("test2.txt");

//string * sequence = readFASTA("CP001872.fna");
//string * sequence2 = readFASTA("CP003913.fna");
//string FASTA1 = "CP001872.fna";
//string FASTA2 = "CP003913.fna";

string * sequence = readFASTA("hg38_chr14_90Mb_100Mb.fa");
string * sequence2 = readFASTA("mm10_chr12_100Mb_110Mb.fa");
string FASTA1 = "hg38_chr14_90Mb_100Mb.fa";
string FASTA2 = "mm10_chr12_100Mb_110Mb.fa";





int cut1 = (*sequence).size();
int cut2 = (*sequence2).size() + cut1 + 1;

vector<vector<char *> > matchStrings(0);

int len = ((*sequence).size() + 2*(*sequence2).size());       

string * reverseSequence = generateReverseComplement(sequence2);


char* allSeq = new char[len + 3];


int *matches = new int[cut1 - 1];
int *histogram = new int[cut1 - 1];
void reportMatches();
void countLongestMatch(int);
int *sortMatrix = new int[len+3];
int *highestLength = new int;  
int countUCEs();


//main function
int main(){
    auto start = chrono::high_resolution_clock::now();
    outfile.open("HOMEWORK1_TIERNAN_KENNEDY.txt");
    outfile << "Assignment: Assignment1" << endl;
    outfile << "Name: Tiernan Kennedy" << endl;
    outfile << "Email: tiernan7@cs.washington.edu"  << endl;
    outfile << "Language: C++" << endl;
    outfile << endl;   
    outfile <<"*="<< (*sequence).size() << endl;
    outfile <<"A=" << countBases(sequence, 'A') << endl;
    outfile << "C="<<countBases(sequence, 'C') << endl;
    outfile << "G="<<countBases(sequence, 'G') << endl;
    outfile << "T="<<countBases(sequence, 'T') << endl;
    outfile << "N="<<countBases(sequence, 'N') << endl;
    outfile << endl;
    outfile <<"*="<< (*sequence2).size() << endl;
    outfile <<"A=" << countBases(sequence2, 'A') << endl;
    outfile << "C="<<countBases(sequence2, 'C') << endl;
    outfile << "G="<<countBases(sequence2, 'G') << endl;
    outfile << "T="<<countBases(sequence2, 'T') << endl;
    outfile << "N="<<countBases(sequence2, 'N') << endl;
    outfile << endl;
    
    //consolodating input strings
    memcpy(allSeq,(*sequence).c_str(),sizeof(char)*((*sequence).length() + 1));
    delete sequence;

    memcpy(&allSeq[cut1 + 1],(*sequence2).c_str(),sizeof(char)*((*sequence2).length() + 1));
    delete sequence2;


    memcpy(&allSeq[cut2 + 1],(*reverseSequence).c_str(),sizeof(char)*((*reverseSequence).length() + 1));
    delete reverseSequence;


    //Generate matrix of pointers to the positions in the sequence matrix
    for (int i = 0; i < len + 3; i++){
        sortMatrix[i] = i;
    }
    
    //quicksort the sequence based of interegers based on the lexographic comparison of the cooresponding sequence array positions
    qsort(sortMatrix,len + 3, sizeof(int), compare);

    


    //Calculated and prints the histogram and other relevant match data
    reportMatches();
    
    

    outfile << endl;
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    outfile << "Runtime: " << duration.count() << " seconds " << endl;
    outfile.close();

    return 0;
}

//Functions called from main that require variables allocated dynamically in main


int compare (const void * a, const void * b){
    int ind1 = *(int*)a;
    int ind2 = *(int*)b;
    const char *str1 = (const char *)&allSeq[ind1];
    const char *str2 = (const char *)&allSeq[ind2];
    return strcmp(str1, str2);

}

int getStrand(int origIndex){
    if (origIndex < cut1){
        return 1;
    }
    else if (origIndex < cut2){
        return 2;
    }
    else if (origIndex == cut1 || origIndex == cut2 || origIndex == len + 3){
        return -1;
    }
    else {
        return 3;
    }
}

int matchLength(char *seq1, char *seq2){
    int i = 0;
    while (seq1[i] == seq2[i]){
        i++;
    }
    return i;
}


string getFasta(int i){
    if (i == 1){
        return FASTA1;
    }
    else {
        return FASTA2;
    }
}

string getStrandText(int i){
    if ( i == 1 ){
        return "Forward";
    }
    else if (i == 2){
        return "Forward";
    }
    else {
        return "Reverse";
    }
} 



void reportMatches(){
    outfile << "Match Length Histogram:" << endl;

    

    for (int i = 0;i < len + 3; i++){
        int origIndex = sortMatrix[i];
        if (getStrand(origIndex) == 1){
            countLongestMatch(i);
         }
    }
    
    
    int histogramEntries = 0;
    for (int i = 0; i < cut1 - 1; i++){
        if (histogram[i] != 0){
            histogramEntries = histogramEntries + histogram[i];
            outfile << i << " " << histogram[i] << endl;
        }
    }
    outfile << "Number of Histogram Entries: " << histogramEntries << endl;
    outfile << endl;
    outfile << "Based on the histogram there are " << countUCEs() << " UCEs in this region." << endl;
    outfile << endl;
    outfile << "Longest match length: " << *highestLength<<endl;
    outfile << "Number of match strings: " << matchStrings.size()  << endl;
    for (int i = 0; i < matchStrings.size(); i++){
        outfile << endl;
        outfile << "Match" << i+1 << ": " <<endl;
        outfile << string(matchStrings[i][0]).substr(0,*highestLength) << endl;
        outfile << "Description:" << endl;
        int position = -1;
        for (int j = 0; j < matchStrings[i].size();j++){
            int index = matchStrings[i][j] - allSeq;
            if (getStrand(index) == 3){
                position = (index - (cut2+1)) - *highestLength + 1;
            }
            else if(getStrand(index) == 2){
                position = index - (cut1 + 1) + 1;
            }
            else {
                position = index + 1;
            }
            outfile << endl;
            outfile << "FATSA: "<< getFasta(getStrand(index)) << endl;
            outfile << "Position: "<< position  <<endl;
            outfile << "Strand :" <<  getStrandText(getStrand(index));
            outfile << endl;

        }
    }
        

    outfile << endl;
    outfile << "Program: " << endl;

}



void countLongestMatch(int sortedPosition){
    int leftMatch = sortedPosition - 1;
    int rightMatch = sortedPosition + 1;
    int *matchPosition = &matches[sortMatrix[sortedPosition]];
    char *match = nullptr;
    while (leftMatch >= 0){
        if (getStrand(sortMatrix[leftMatch]) != 1 && getStrand(sortMatrix[leftMatch]) != -1){
            *matchPosition = matchLength(&allSeq[sortMatrix[leftMatch]], &allSeq[sortMatrix[sortedPosition]]);
            match = &allSeq[sortMatrix[leftMatch]]; 
            break;
        }
        leftMatch--;
    }   
    while (rightMatch < len + 3){
        if (getStrand(sortMatrix[rightMatch]) != 1 && getStrand(sortMatrix[rightMatch]) != -1){
            int secondMatch = matchLength(&allSeq[sortMatrix[rightMatch]], &allSeq[sortMatrix[sortedPosition]]);
            if (secondMatch > *matchPosition){
                *matchPosition = secondMatch;
                match = &allSeq[sortMatrix[rightMatch]]; 
            }     
            break;
        }
        rightMatch++;
    }
    histogram[*matchPosition]++;
    if (*matchPosition == *highestLength){
        if ((string(&allSeq[sortMatrix[sortedPosition]]).substr(0,*highestLength).compare(string(match).substr(0,*highestLength))) != 0){
            *highestLength = *matchPosition;
            vector <char*> mstring(0);
            mstring.push_back(&allSeq[sortMatrix[sortedPosition]]);
            mstring.push_back(match);
            matchStrings.push_back(mstring);
        }
        else{
            bool true1 = false;
            bool true2 = false;
            vector<char *> question = matchStrings.back();
            for (int i = 0; i < question.size() ;i++){
                if (question[i] == match){
                    true1 = true;
                }
                if (question[i] == &allSeq[sortMatrix[sortedPosition]]){
                    true2 = true;
                }
            }
            if (not (true1)){
                    (matchStrings.back()).push_back(match);

            }
            if (not (true2)){
                    (matchStrings.back()).push_back(&allSeq[sortMatrix[sortedPosition]]);
            }
        }
    }
    else if (*matchPosition > *highestLength){
        *highestLength = *matchPosition;
        matchStrings.clear();
        vector<char *> mstring(0);
        mstring.push_back(&allSeq[sortMatrix[sortedPosition]]);
        mstring.push_back(match);
        matchStrings.push_back(mstring);
    }
}    


int countUCEs(){
    int count = 0;
    for (int i = 200; i < cut1-1;i++ ){
        count = count + histogram[i]; 
    }
    return count;
}

