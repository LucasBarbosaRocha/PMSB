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
    cout << "Qtd. Kmers: " << h.getQtdKmers() << endl;
    h.displayHash();
    auto grafo2 = h.dbgToSequenceGraph_2();
    h.displayHash();

    auto m_grafo2 = m.buildMultilayerGraph(grafo2.second, sequence);
    cout << "Qtd. Kmers + especiais: " << h.getQtdKmers() << endl;
    cout << "Qtd. vértices sequence graph: " << grafo2.second.getV() << endl;
    cout << "Qtd. vértices multicamada: " << m_grafo2.getV() << endl;
    m.shortestPath(m_grafo2, m.getInitialNode(), m.getEndNode(), 1);

    /*retorno = m.dijkstra(m_grafo2,  m.getInitialNode(), m.getEndNode());
    auto saida2 = m.mostraMapeamento_2(retorno.first, grafo2.first, h.getKmerSpecialAndKmer(), grafo2.second);
    cout << "Percurso: ";
    for (auto it : saida2.first)
    {
        cout << it << " ";
    }
    cout << "Mapeamento: " << saida2.second << " cost: " << retorno.second << endl;  */

return 0;
}
