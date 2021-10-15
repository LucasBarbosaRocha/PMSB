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
    h.displayHash();

    cout << "Qtd. Kmers: " << h.getQtdKmers() << endl;
    cout << "Qtd. vértices sequence graph: " << grafo.second.getV() << endl;
    cout << "Qtd. vértices multicamada: " << m_grafo.getV() << endl;
    m.shortestPath(m_grafo, m.getInitialNode(), m.getEndNode(), 1);
    //auto retorno = m.dijkstra(m_grafo, m.getInitialNode(), m.getEndNode());
    //auto saida = m.mostraMapeamento(retorno.first, grafo.first, grafo.second);   
    //cout << "Percurso: ";
    //for (auto it : saida.first)
    //{
    //    cout << it << " ";
    //}
    //cout << "Mapeamento: " << saida.second << " cost: " << retorno.second << endl;  

return 0;
}
