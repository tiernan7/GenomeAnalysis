#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;


class vertex(string label, int startOrEnd){
    assert(label.size() <= 10);
    assert(startOrEnd == 0 or startOrEnd == 1 or startOrEnd == 2);
    this.label = label;
    this.startOrEnd = startOrEnd;
};


def readFASTA()
