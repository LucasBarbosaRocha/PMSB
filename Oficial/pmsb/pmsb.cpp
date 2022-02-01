#include <iostream>
//#include "../utils/myDBGgraph.cpp"
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"
#include <queue>
#include <bits/stdc++.h>

int main(int argc, char *argv[])
{
    //string nomeArquivo = "../archives/kmers.txt";
    Marschall m;  
    MyUtils utils;
    if (utils.verifyData(argc, argv) == 1)
        exit(0);

    Hash h(utils.k); 
    utils.readSequence(utils.nameSequenceArchive);
    h.populateGraph(utils.nameArchive, false); 
    //h.displayHash();

    if(utils.typeGraph == 0)
    {
        SequenceGraph traditionalGraph = h.dbgToTraditionalSequenceGraph();
        //traditionalGraph.printGraph();
        SequenceGraph m_traditionalGraph = m.buildMultilayerGraph(traditionalGraph, utils.sequence);
        auto retorno = m.dijkstra(m_traditionalGraph, m.getInitialNode(), m.getEndNode());
        auto saida = m.showTraditionalMapping(retorno.first, h, traditionalGraph);
        cout << "Traditional." << endl;
        cout << "Cost: " << retorno.second << endl; 
        cout << "Mapping:" << saida.second << endl;    
    } else {
        SequenceGraph simplifedGraph = h.dbgToSimplifiedSequenceGraph();
        //simplifedGraph.printGraph();
        SequenceGraph m_simplifedGraph = m.buildMultilayerGraph(simplifedGraph, utils.sequence);
        auto retorno = m.dijkstra(m_simplifedGraph, m.getInitialNode(), m.getEndNode());
        auto saida = m.showSimplifiedMapping(retorno.first, h, simplifedGraph);
        cout << "Simplified." << endl;
        cout << "Cost: " << retorno.second << endl; 
        cout << "Mapping:" << saida.second << endl;    
    }


return 0;
}