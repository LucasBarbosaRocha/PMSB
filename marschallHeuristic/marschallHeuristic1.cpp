#include <iostream>
#include "marschall.cpp"
#include <queue>
#include <bits/stdc++.h>

int main (void)
{
    string sequence = "AAGACGTGTAGACGTGA", nomeArquivo = "kmers4.txt";
    // AAGACGTGTAGACGTGA
    map<string, int> kmersInTheGraph;
    vector<int> posicoesValidas;
    int k = 3;
    //cin >> sequence;
    //cin >> nomeArquivo;
    //cin >> k;
    Hash h(k);   
    Marschall m;
    h.populateGraph(nomeArquivo, false);

    for (int i = 0; i < sequence.length() - (k - 1); i++)
    {
        if (h.contains(sequence.substr(i,k)))
        {
            cout << sequence.substr(i,k) << endl;
            posicoesValidas.push_back(i);
        }
    }

    auto grafoSequencias = h.dbgToSequenceGraph_1();
    
    // mapeamento 
    cout << sequence << endl;
    if (posicoesValidas.size() > 0)
    {
        if (posicoesValidas[0] != 0)
        {
            cout << "(-," << posicoesValidas[0] << ") ";
        }

        for(int i = 0; i < posicoesValidas.size() - 1; i++)
        {     
            cout << "(" << posicoesValidas[i] << "," << posicoesValidas[i + 1] << ") "; 
            int dif = posicoesValidas[i + 1] -  posicoesValidas[i]; 
            if (dif > 1)
            {
                string sequence_aux = sequence.substr(posicoesValidas[i], k + dif);
                cout << "vamos procurar " << sequence_aux << endl;
                auto grafo_multicamada = m.buildMultilayerGraph(grafoSequencias.second, sequence_aux);
                auto retorno = m.dijkstra(grafo_multicamada, m.getInitialNode(), m.getEndNode());
                auto mapeamento = m.mostraMapeamento(retorno.first, grafoSequencias.first, grafoSequencias.second);    
                cout << mapeamento << " com custo " << retorno.second << endl;            
            }
        }

        if (posicoesValidas[posicoesValidas.size() - 1] != sequence.length() - k)
        {
            cout << "(" << posicoesValidas[posicoesValidas.size() - 1] << ",-)" << endl;
        }
    } else
    {
        cout << "Nenhum kmer encontrado no grafo. A sequência vai ser toda mapeada." << endl;
        auto grafo_multicamada = m.buildMultilayerGraph(grafoSequencias.second, sequence);
        auto retorno = m.dijkstra(grafo_multicamada, m.getInitialNode(), m.getEndNode());
        auto mapeamento = m.mostraMapeamento(retorno.first, grafoSequencias.first, grafoSequencias.second);    
        cout << mapeamento << " com custo " << retorno.second << endl;            
    }


    return 0;
}