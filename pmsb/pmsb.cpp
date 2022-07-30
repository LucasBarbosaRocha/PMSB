#include <iostream>
#include <fstream>
//#include "../utils/myDBGgraph.cpp"
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"
#include <queue>
#include <bits/stdc++.h>

int main(int argc, char *argv[])
{
    //string nomeArquivo = "../archives/kmers.txt";
    string line;
    Marschall m;  
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
        if(utils.typeGraph == 0)
        {
            m.buildMultilayerGraph(h.sequenceGraph, utils.sequence);
            auto retorno = m.dijkstra(m.m_sequenceGraph, m.getInitialNode(), m.getEndNode());
            auto saida = m.showTraditionalMapping(retorno.first, h, h.sequenceGraph);
            //cout << "Traditional." << endl;
            cout << saida.second << endl;    
            cout << "Cost: " << retorno.second << endl; 
            
        } else {
            m.buildMultilayerGraph(h.sequenceGraph, utils.sequence);
            auto retorno = m.dijkstra(m.m_sequenceGraph, m.getInitialNode(), m.getEndNode());
            auto saida = m.showSimplifiedMapping(retorno.first, h, h.sequenceGraph);
            //cout << "Simplified." << endl;
            cout << saida.second << endl;    
            cout << "Cost: " << retorno.second << endl; 
            m.m_sequenceGraph.deleteGraph();
        }
    }

return 0;
}


