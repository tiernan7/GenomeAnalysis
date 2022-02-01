#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "main.h"
#include <regex>
#include <stdexcept>
#include <boost/tokenizer.hpp>
#include <sstream>

using namespace std;

//NOTES: This program assumes that the input graph is topologically sorted, and is thus acyclic


//A class representing a node/vertex in a directed graph
Vertex::Vertex(int lab, int sOrE){
    assert(to_string(lab).length() <= 10);
    assert(sOrE == 0 or sOrE == 1 or sOrE == -1);
    label = lab;
    startOrEnd = sOrE;
}
    
bool Vertex::isStart(){
    if (startOrEnd == 1){
        return true;
    } 
    else{
        return false;
    }
}
   
bool Vertex::isEnd(){
    if (startOrEnd == -1){
        return true;
    }
    else{
        return false;
    }
}

int Vertex::getNext(){
      return next;
}

void Vertex::setNext(int label){
        next = label;
}

int Vertex::getLabel(){
        return label;
}

void Vertex::addParent(int label){
    parents.push_back(label);
}
    
bool Graph::multipleStarts(){
    int numberOfStarts = 0;
    for (int i = 0; i < nodeList->size(); i++){
        if (nodeList->at(i)->isStart()){
            numberOfStarts++;
        }
    }
    return (numberOfStarts>1);
}
    
bool Graph::multipleEnds(){
    int numberOfEnds = 0;
    for (int i = 0; i < nodeList->size(); i++){
        if (nodeList->at(i)->isEnd()){
            numberOfEnds++;
        }
    return (numberOfEnds > 1);
    }
}

bool Graph::isValidGraph(){
    bool returnValue = false;
    returnValue = not(multipleStarts() or multipleEnds());
    return returnValue;
}

void Graph::processLine(string line){
    
    boost::char_separator<char> sep(" ");
    boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
    vector<string> words;

    for (const auto &t: tokens){
        words.push_back(t);
    } 
    

    int startOrEnd = 0;
    if (words[0] == "V"){

        stringstream l(words[1]);
        int label;
        l >> label;
        for (auto &t : words){
            if (t == "Start"){
                startOrEnd = 1;
            }        
            if (t == "End"){
            startOrEnd = -1;
           } 
        }
        


        Vertex * v = new Vertex(label, startOrEnd);
        
        nodeList->push_back(v);
        }
        
    else if (words[0] == "E"){
        stringstream start(words[2]);
        stringstream end(words[3]);
        int s = 0;
        int e = 0;
        start >> s;
        end >> e;
        (*nodeList)[int(s)]->setNext(s);
        (*nodeList)[int(e)]->addParent(s);    
     }

    else {
        throw invalid_argument( "File does not exist" );
    }
}

int Graph::countNodes(){
    return nodeList->size();
}

Graph::Graph(string fileName){
    ifstream graphFile;
    graphFile.open(fileName);

    if (graphFile.is_open()){
        string line;
        while (getline(graphFile,line)){
            processLine(line);
        }
    }
    else{
        throw invalid_argument( "File does not exist" );
    }
    assert(isValidGraph());
}


int main(){
    Graph * input = new Graph("Graph.txt");
    cout << input->countNodes() << endl;
    return 0;
}
