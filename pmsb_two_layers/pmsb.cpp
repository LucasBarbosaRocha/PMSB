#include <iostream>
#include <fstream>
//#include "../utils/myDBGgraph.cpp"
#include "../utils/marschallTwoLayers.cpp"
#include "../utils/myUtils.cpp"
#include <queue>
#include <bits/stdc++.h>


int main(int argc, char *argv[])
{
    //string nomeArquivo = "../archives/kmers.txt";
    string line;
    MarschallTwoLayers m;  
    MyUtils utils;
    if (utils.verifyData(argc, argv) == 1)
        exit(0);

    Hash h(utils.k); 
    ifstream file(utils.nameSequenceArchive);
    // utils.readSequence(utils.nameSequenceArchive);
    h.populateGraph(utils.nameArchive, false); 
    
   if(utils.typeGraph == 0)
        h.dbgToSimplifiedSequenceGraph(0);
    else
        h.dbgToSimplifiedSequenceGraph(1);

    while(getline(file, line))
    {
        //utils.readSequence(utils.nameSequenceArchive);  
        getline(file, line);
        cout << "Size L.Read " << line.size() << endl;
        utils.sequence = line;
        // mapeamento
        int cost = m.toCalculateTheCost(h.sequenceGraph, utils.sequence);
        if (cost != INT_MAX)
            cout << "Cost: " << cost << endl; 
        else {
            cout << "Sem mapeamento\nCost: " << cost << endl;
        }
    }

return 0;
}


