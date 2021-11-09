#include <iostream>
#include "marschall.cpp"
#include <queue>
#include <bits/stdc++.h>

string mapeamento(Hash h, pair<unordered_map<int, string>,SequenceGraph> grafoSequencias, string sequence, int k)
{
    cout << sequence << endl;
    int posicao = 0;
    string resposta = "";
    vector<int> posicoesValidas;
    Marschall m;

    for (int i = 0; i < sequence.length() - (k - 1); i++)
    {
        if (h.contains(sequence.substr(i,k)))
        {
            posicoesValidas.push_back(i);
        }
    }

    if (posicoesValidas.size() > 0)
    {
        if (posicoesValidas[0] != 0)
        {
            //cout << "(-," << posicoesValidas[0] << ") ";
            string sequence_aux = sequence.substr(0, k + posicoesValidas[0]);
            auto grafo_multicamada = m.buildMultilayerGraph(grafoSequencias.second, sequence_aux);
            auto retorno = m.dijkstra(grafo_multicamada, m.getInitialNode(), m.getEndNode());
            auto mapeamento = m.mostraMapeamento(retorno.first, grafoSequencias.first, grafoSequencias.second);    
            resposta = mapeamento.second.substr(0,posicoesValidas[0]) + sequence.substr(posicoesValidas[0], k); 
            posicao++;
        }

        for(int i = 0; i < posicoesValidas.size() - 1; i++)
        {     
            //cout << "(" << posicoesValidas[i] << "," << posicoesValidas[i + 1] << ") "; 
            int dif = posicoesValidas[i + 1] -  posicoesValidas[i]; 
            if (dif > 1)
            {
                string sequence_aux = sequence.substr(posicoesValidas[i], k + dif);
                auto grafo_multicamada = m.buildMultilayerGraph(grafoSequencias.second, sequence_aux);
                auto retorno = m.dijkstra(grafo_multicamada, m.getInitialNode(), m.getEndNode());
                auto mapeamento = m.mostraMapeamento(retorno.first, grafoSequencias.first, grafoSequencias.second);    
                if (mapeamento.second.length() >= k)      
                {
                    resposta = resposta + mapeamento.second.substr(k, mapeamento.second.length() - (2 * k));   
                    //resposta = resposta + mapeamento.substr(mapeamento.length()- (k), k);   
                    resposta = resposta + sequence.substr(posicoesValidas[i+1], k);                 
                    posicao = 1;
                }
            }else{
                if (posicao == 0)
                    resposta = resposta + sequence.substr(posicoesValidas[i], k);
                else
                    resposta = resposta + sequence.substr(posicoesValidas[i] + (k), 1);
                posicao++;
            }
        }

        if (posicoesValidas[posicoesValidas.size() - 1] != sequence.length() - k)
        {
            //cout << "(" << posicoesValidas[posicoesValidas.size() - 1] << ",-)" << endl;
            int dif = sequence.length() - posicoesValidas[posicoesValidas.size() - 1];
            if (dif > 1)
            {    
                string sequence_aux = sequence.substr(posicoesValidas[posicoesValidas.size() - 1], k + dif);
                auto grafo_multicamada = m.buildMultilayerGraph(grafoSequencias.second, sequence_aux);
                auto retorno = m.dijkstra(grafo_multicamada, m.getInitialNode(), m.getEndNode());
                auto mapeamento = m.mostraMapeamento(retorno.first, grafoSequencias.first, grafoSequencias.second);    
                if (mapeamento.second.length() >= k)      
                    resposta = resposta + mapeamento.second.substr(k, mapeamento.second.length() - k);  
            }
        }
    } else
    {
        cout << "Nenhum kmer encontrado no grafo. A sequência vai ser toda mapeada." << endl;
        auto grafo_multicamada = m.buildMultilayerGraph(grafoSequencias.second, sequence);
        auto retorno = m.dijkstra(grafo_multicamada, m.getInitialNode(), m.getEndNode());
        auto resposta = m.mostraMapeamento(retorno.first, grafoSequencias.first, grafoSequencias.second);    
        return resposta.second;
    }
    return resposta;
}

int main (void)
{
    string sequence = "TGCAAAACTATAAAAATTTTGGAACGAGTCAAAGGGTAGGTCTTGGGAATACTACTTCATCCGGTTCTTG", nomeArquivo = "/home/lucas/Documentos/sequencias/sequencias100.fasta";
    // AAGACGTGTAGACGTGA
    map<string, int> kmersInTheGraph;
    int k = 10;
    Hash h(k);   
    h.populateGraph(nomeArquivo, false);
    auto grafoSequencias = h.dbgToSequenceGraph_1();    
    // mapeamento
    auto retorno = mapeamento(h, grafoSequencias, sequence, k); 
    cout << "Mapeamento " << retorno << endl;
    return 0;
}