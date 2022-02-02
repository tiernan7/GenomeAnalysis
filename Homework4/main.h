#ifndef MAIN
#define MAIN

//A struct representing an edge
//Each edge is stored to be stored within a vertex object
struct Edge{
    int otherEnd;
    double weight;
};

//A function for logging to a File
//


//A class representing a vertex on a graph
class Vertex{
    private:
        int label;
        int startOrEnd;
        std::vector<Edge> children;
        std::vector<Edge> parents;
        int highestParent;
        double highestWeight;

    public:
        Vertex(int label, int sOrE);
        bool isStart();
        bool isEnd();
        void addChild(int label, double weight);
        int getLabel();
        void addParent(int label, double weight);
        std::vector<Edge> getParents();
        int getHighestParent();
        void setHighestParent(int highest);
        double getHighestWeight();
        void setHighestWeight(double weight);
};

//A class representing an entire graph with member functions for graph search
//  and manipulation
class Graph{
    public:
        std::vector<Vertex*> * nodeList = new std::vector<Vertex*>;
        int startNode = -999;
        int endNode = -999;
    
        bool multipleStarts();
        bool multipleEnds();
        bool isValidGraph();
        void processLine(std::string line);
        int countNodes();
        
        Graph(std::string fileName);

        std::vector<int> * getHighestPath(int start, int end);

};


#endif

