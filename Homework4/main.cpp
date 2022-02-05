#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "main.h"
#include <regex>
#include <stdexcept>
#include <boost/tokenizer.hpp>
#include <sstream>
#include <chrono>
#include <cmath>

using namespace std;

//NOTES: This program assumes that the input graph is topologically sorted, and is thus acyclic


//A class representing a node/vertex in a directed graph
Vertex::Vertex(int lab, int sOrE){
    assert(to_string(lab).length() <= 10);
    assert(sOrE == 0 or sOrE == 1 or sOrE == -1);
    label = lab; //An integer label for the vertex
    startOrEnd = sOrE; //-1 when an end, 0 when neither, 1 when a start
}
    
bool Vertex::isStart(){ //Returns true if the vertex is a start vertex
    if (startOrEnd == 1){
        return true;
    } 
    else{
        return false;
    }
}
   
bool Vertex::isEnd(){//Return true if the vertex is an end vertex
    if (startOrEnd == -1){
        return true;
    }
    else{
        return false;
    }
}

vector<Edge> Vertex::getChildren(){
    return children;
}

void Vertex::addChild(int label,double weight, string l){//Adds an edge directed at a child with a specified weight
        struct Edge newEdge;
        newEdge.otherEnd = label;
        newEdge.weight = weight;
        newEdge.label = label;
        children.push_back(newEdge);
}

int Vertex::getLabel(){ //Returns the label of the vertex as an int
        return label;
}

void Vertex::addParent(int label, double weight, string l){//Adds an edge directed from a parent to the current node with a specified weight
    struct Edge newEdge;
    newEdge.otherEnd = label;
    newEdge.weight = weight;
    newEdge.label = l;
    parents.push_back(newEdge);
}

vector<Edge> Vertex::getParents(){//Returns the vector of incoming edges
    return parents;
}
    
int Vertex::getHighestParent(){//Returns the label of the incoming path with highest weight
    return highestParent;
}

void Vertex::setHighestParent(int highest){//Sets the value of the incoming path with the highest weight (or the current node)
    highestParent = highest;
}

double Vertex::getHighestWeight(){//Returns the value of the weight of the highest incoming edge, or 0
    return highestWeight;
}

void Vertex::setHighestWeight(double weight){//Sets the value of the highest weight
    highestWeight = weight;
}

void Vertex::setHighestChild(int c){
    highestChild = c;
}

bool Graph::multipleStarts(){//Returns true if multiple start sites exist in the graph
    int numberOfStarts = 0;
    for (int i = 0; i < nodeList->size(); i++){
        if (nodeList->at(i)->isStart()){
            numberOfStarts++;
        }
    }
    return (numberOfStarts>1);
}


    
bool Graph::multipleEnds(){//Returns true if multiple end sits exist in the graph
    int numberOfEnds = 0;
    for (int i = 0; i < nodeList->size(); i++){
        if (nodeList->at(i)->isEnd()){
            numberOfEnds++;
        }
    }
    
    return (numberOfEnds > 1);
}

bool Graph::isValidGraph(){//Checks to make sure there are not multiple starts or ends
    bool returnValue = false;
    returnValue = not(multipleStarts() or multipleEnds());
    return returnValue;
}

void Graph::processLine(string line){//Creates a node or edge from a given line of a graph input file
    
    boost::char_separator<char> sep(" "); //tokenizer to break each line into words
    boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
    vector<string> words;

    for (const auto &t: tokens){
        words.push_back(t);
    } 

    int startOrEnd = 0;
    if (words[0] == "V"){//Handles the case where the line is a vertex
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
        //creates the vertex and adds it to the nodes list
        nodeList->push_back(v);
        }
        
    else if (words[0] == "E"){//Handles the case where the line is an edge
        stringstream start(words[2]);
        stringstream end(words[3]);
        stringstream w(words[4]);
        string label = words[1];
        double wt;
        w >> wt;
        int s = 0;
        int e = 0;
        start >> s;
        end >> e;
        (*nodeList)[int(s)]->addChild(s,wt,label);
        (*nodeList)[int(e)]->addParent(s,wt,label);
        
        
        
     }

    else {
        throw invalid_argument( "File does not exist" );
    }
}

int Graph::countNodes(){ //Returns the number of nodes 
    return nodeList->size();
}

Graph::Graph(string fileName){//Initializes a graph from a file
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


//Function to set start node of graph, apply only after isValidGraph() has been berified
void Graph::setStart(){
    for (auto &n : *nodeList){
        if (n->isStart()){
            startNode = n->getLabel();
        }
    }
}

//Function to set end node of graph, apply only after isValidGraph() has been berified
void Graph::setEnd(){
    for (auto &n : *nodeList){
        if (n->isEnd()){
            endNode = n->getLabel();
        }
    }
}


void Graph::setPathStart(int s){
    pathStart = s;
}

void Graph::setPathEnd(int e){
    pathEnd = e;
}

void Graph::setPath(string p){
    pth = p;
}

void Graph::setScore(double s){
    score = std::ceil(100.0*s)/100.0;
}


//A function that returns the highest weight path on a graph as a vector of labels
vector<int> * Graph::getHighestPath(int start, int end){
    vector<int> * returnPath = new vector<int>; //the longest path to be returned
    if (start == -999 and end == -999){ //case with no constraints

        for (auto &v : *nodeList){ //for each node in the graph

            if (v->getParents().size() == 0){ //if depth i s0...

                v->setHighestParent(-1); //highest parent is undefined
                v->setHighestWeight(0.0); //highest weight must be 0

            }
            else{ //else we must have visited the parents, so...
                double highestWeight = -numeric_limits<double>::infinity(); //initialize for finding max
                double highestLabel = -1;

                for (auto &p : v->getParents()){ //for each parent of the node
                    if (p.weight > highestWeight){ //find highest weight parent
                        highestWeight = p.weight;
                        highestLabel = p.otherEnd;
                    }

                }
                if (highestWeight + nodeList->at(highestLabel)->getHighestWeight()  < 0){
                    highestWeight = 0.0;
                    highestLabel = -1;
                    v->setHighestParent(highestLabel); 
                    v->setHighestWeight(highestWeight);
                }
                else{
                    v->setHighestParent(highestLabel); 
                    v->setHighestWeight(nodeList->at(highestLabel)->getHighestWeight() + highestWeight);
                }
            }
            
        }
        string s = "";
        double highest = -numeric_limits<double>::infinity(); //varibles to keep track of the longest path
        int highestLabel = -9999;

         for (auto &v : *nodeList){
            if (v->getHighestWeight() > highest){
                highest = v->getHighestWeight();
                highestLabel = v->getLabel();
         }
        }

        int pos = highestLabel;
        returnPath->push_back(pos);
        do{
            for (auto &p : nodeList->at(pos)->getParents()){
                if (p.otherEnd == nodeList->at(pos)->getHighestParent()){
                    s += p.label;
                }
            }

            pos = nodeList->at(pos)->getHighestParent();
            returnPath->push_back(pos);
        
        } while ((*nodeList).at(pos)->getHighestParent() != -1);


        reverse(returnPath->begin(), returnPath->end());
        reverse(s.begin(),s.end());

        setScore(highest);
        setPathStart(*returnPath->begin());

        setPathEnd(returnPath->back());
        setPath(s);

        return returnPath;
    }

    else{ //Case where start and end nodes are set
        
        //

        setPathStart(startNode);
        setPathEnd(endNode);
        string p = ""; 
        for (int i = startNode; i <= endNode; i++){
            if (nodeList->at(i)->getParents().size() == 0){
                if (i == startNode){
                    nodeList->at(i)->setHighestWeight(0);
                    nodeList->at(i)->setHighestParent(0);
                }
                else {
                    nodeList->at(i)->setHighestWeight(0.0);
                    nodeList->at(i)->setHighestParent(-99);

                }
                
            }
            else{
                double highestLabel = -1;
                double highestWeight = -numeric_limits<double>::infinity(); //initialize for finding max
                
               for (auto &n : nodeList->at(i)->getParents()){
                    cout<< "parent label: " << n.label << endl;
                   if (nodeList->at(n.otherEnd)->getHighestParent() == -99){
                        cout << "skipped" << highestWeight <<endl;
                       cout<< n.weight + nodeList->at(n.otherEnd)->getHighestWeight() <<endl;
                       cout <<".."<<endl;
                       cout <<endl;
                       ;
                   }
                   else if (n.weight + nodeList->at(n.otherEnd)->getHighestWeight() > highestWeight){
                       cout << highestWeight << endl;
                       cout << nodeList->at(n.otherEnd)->getHighestWeight() <<endl;
                       cout << n.weight << endl;
                       cout<< n.weight + nodeList->at(n.otherEnd)->getHighestWeight() <<endl;          highestWeight = n.weight + nodeList->at(n.otherEnd)->getHighestWeight();          highestLabel = n.otherEnd;

                    cout <<".."<<endl;
                

                        
                  }
                   else{
                       cout << highestWeight<<endl;
                       cout << nodeList->at(n.otherEnd)->getHighestWeight() <<endl;
                       cout << n.weight << endl;
                       cout<< "too low "<< n.weight + nodeList->at(n.otherEnd)->getHighestWeight() <<endl;
                   }
                        
                }
                cout << "final: " << highestWeight << endl;
                if (nodeList->at(i)->getParents().size() != 0){
                    cout << "gh" << endl;
                    nodeList->at(i)->setHighestWeight(highestWeight);
                    nodeList->at(i)->setHighestParent(highestLabel);
                }
                cout << "...."<<endl;
                cout<<endl;    

              }
        }
        setScore(nodeList->at(endNode)->getHighestWeight());
    }
}




string Graph::getScore(){
    string returnString;
    returnString = to_string(std::ceil(100.0)*score/100); 
   
    return returnString;
}


string Graph::getPathStart(){
    string returnString;
    returnString = to_string(pathStart);
    return returnString;
}

string Graph::getPathEnd(){
    string returnString;
    returnString = to_string(pathEnd);
    return returnString;
}

string Graph::getPath(){
    return pth;
}


string Graph::summarize(){
    string returnString = "";
    returnString += "Score: " + getScore() + "\n";
    returnString += "Begin: " + getPathStart() + "\n";
    returnString += "End: " + getPathEnd() + "\n";
    returnString += "Path: " + getPath() + "\n";
    returnString += "\n";

    return returnString;

}



int main(){
    cout.precision(2);
    Graph * part1 = new Graph("graph.txt");
    part1->getHighestPath(part1->startNode, part1->endNode);
    cout << part1->summarize();
   
    part1->setStart();
    part1->setEnd();
    part1->getHighestPath(part1->startNode, part1->endNode);
    cout << part1->summarize();

    Graph * part3 = new Graph("graphGCF_000967895.1_ASM96789v1_genomic.fna");
    part3->getHighestPath(part3->startNode, part3->endNode);
    cout << part3->summarize();

    
    return 0;
}
