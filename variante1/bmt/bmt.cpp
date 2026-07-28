#include <iostream>
#include <fstream>
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"
#include <queue>
#include <bits/stdc++.h>

constexpr int kUseHeuristic = 1;

int main(int argc, char *argv[])
{
    string line;
    Marschall m;
    MyUtils utils;
    if (utils.verifyData(argc, argv) == 1)
        exit(0);

    Hash h(utils.k);
    ifstream file(utils.nameSequenceArchive);
    h.populateGraph(utils.nameArchive, false);

    if (utils.typeGraph == 0)
        h.dbgToTraditionalSequenceGraph(0, kUseHeuristic);
    else
        h.dbgToSimplifiedSequenceGraph(0);

    while (getline(file, line))
    {
        getline(file, line);
        cout << "Size L.Read " << line.size() << endl;
        utils.sequence = line;

        m.buildMultilayerGraph(h.sequenceGraph, utils.sequence);
        auto retorno = m.dijkstra(m.m_sequenceGraph, m.getInitialNode(), m.getEndNode());

        if (utils.typeGraph == 0)
        {
            auto saida = m.showTraditionalMapping(retorno.first, h, h.sequenceGraph);
            cout << saida.second << endl;
        }
        else
        {
            auto saida = m.showSimplifiedMapping(retorno.first, h, h.sequenceGraph);
            cout << saida.second << endl;
            m.m_sequenceGraph.deleteGraph();
        }
        cout << "Cost: " << retorno.second << endl;
    }

    return 0;
}
