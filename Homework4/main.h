#ifndef MAIN
#define MAIN



//A struct representing an edge
//Each edge is stored to be stored within a vertex object
struct Edge{
    std::string label;
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
        int highestParent = -99;
        double highestWeight;
        int highestChild;

    public:
        Vertex(int label, int sOrE);
        bool isStart();
        bool isEnd();
        void addChild(int label, double weight, std::string l);
        int getLabel();
        void addParent(int label, double weight, std::string l);
        std::vector<Edge> getParents();
        int getHighestParent();
        void setHighestParent(int highest);
        double getHighestWeight();
        void setHighestWeight(double weight);
        std::vector<Edge> getChildren();
        void setHighestChild(int c);
};

//A class representing an entire graph with member functions for graph search
//  and manipulation
class Graph{
    private:
        int pathStart;
        int pathEnd;
        double score;
        std::string pth = "";
    public:
        std::vector<Vertex*> * nodeList = new std::vector<Vertex*>;
        int startNode = -999;
        int endNode = -999;
    
        bool multipleStarts();
        bool multipleEnds();
        bool isValidGraph();
        void processLine(std::string line);
        int countNodes();
        void setStart();
        void setEnd();
        std::string summarize();
        std::string getScore();
        std::string getPathStart();
        std::string getPathEnd();
        std::string getPath();
        void setScore(double s);
        void setPathStart(int s);
        void setPathEnd(int e);
        void setPath(std::string p);

        
        Graph(std::string fileName);

        std::vector<int> * getHighestPath(int start, int end);

};


#endif

