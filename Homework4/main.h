#ifndef MAIN
#define MAIN

class Vertex{
    private:
        int label;
        int startOrEnd;
        int next;
        std::vector<int> parents;

    public:
        Vertex(int label, int sOrE);
        bool isStart();
        bool isEnd();
        int getNext();
        void setNext(int label);
        int getLabel();
        void addParent(int label);
};

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
};


#endif

