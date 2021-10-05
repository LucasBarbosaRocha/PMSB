#include <iostream>
#include "marschall.cpp"
#include <queue>
#include <bits/stdc++.h>

int main(int argc, char *argv[])
{
    /* string sequence = "CGA";
    int k = 3;
    string nomeArquivo = "kmers4.txt"; */

    if(verificaEntrada(argc, argv) == 1)
        exit (0);
    
    string aux = "";
        
    /* insert the kmers into the hash table */
    Hash h(k);   
    Marschall m;  
    h.populateGraph(nameArchive, false); 
    auto grafo = h.dbgToSequenceGraph_1();
    auto m_grafo = m.buildMultilayerGraph(grafo.second, sequence);
    auto retorno = m.dijkstra(m_grafo, m.getInitialNode(), m.getEndNode());
    auto saida = m.mostraMapeamento(retorno.first, grafo.first, grafo.second);   
    cout << "Mapeamento: " << saida << " cost " << retorno.second << endl; 

    auto grafo2 = h.dbgToSequenceGraph_2();
    auto m_grafo2 = m.buildMultilayerGraph(grafo2.second, sequence);
    aux = "";
    retorno = m.dijkstra(m_grafo2,  m.getInitialNode(), m.getEndNode());
    auto saida2 = m.mostraMapeamento(retorno.first, grafo2.first, grafo2.second);
    cout << "Mapeamento: " << saida2 << " cost: " << retorno.second << endl;  

return 0;
}
