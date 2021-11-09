#include <iostream>
#include "marschall.cpp"
#include <queue>
#include <bits/stdc++.h>

string mapeamento(Hash h, string sequence, int k)
{
    cout << sequence << endl;
    int posicao = 0, length, aux;
    float x = 0.5;
    string resposta = "";
    vector<int> posicoesValidas;
    list<string> kmer_lista_aux;
    Marschall m;
    Hash dbg_gap(k);

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
            length = k * x, aux = 0;

            string kmer_cauda = sequence.substr(posicoesValidas[0], k);

            dbg_gap.insertKmer(kmer_cauda);   
            // buscando os kmers do gap

            for(int j = 0; j < posicoesValidas[0]; j++)
            {

                string kmer_aux = sequence.substr(j, k);
                kmer_lista_aux = h.compareKmersWithGraph(kmer_aux, x * k); 
                for (auto aux : kmer_lista_aux)
                    dbg_gap.insertKmer(aux);
            }  

            auto sequence_graph_aux = dbg_gap.dbgToSequenceGraph_1();    
        
            string sequence_aux = sequence.substr(0, k + posicoesValidas[0]);
            auto grafo_multicamada = m.buildMultilayerGraph(sequence_graph_aux.second, sequence_aux);
            auto retorno = m.dijkstra(grafo_multicamada, m.getInitialNode(), m.getEndNode());
            auto mapeamento = m.mostraMapeamento(retorno.first, sequence_graph_aux.first, sequence_graph_aux.second);    
            resposta = mapeamento.second.substr(0,posicoesValidas[0]) + sequence.substr(posicoesValidas[0], k); 
            posicao++;
        }

        for(int i = 0; i < posicoesValidas.size() - 1; i++)
        {     
            //cout << "(" << posicoesValidas[i] << "," << posicoesValidas[i + 1] << ") "; 
            int dif = posicoesValidas[i + 1] -  posicoesValidas[i]; 
            if (dif > 1)
            {
                length = dif * x, aux = 0;
                string kmer_cabeca = sequence.substr(posicoesValidas[i], k);
                string kmer_cauda = sequence.substr(posicoesValidas[i+1], k);
                dbg_gap.deleteHash();
                dbg_gap.insertKmer(kmer_cabeca);
                dbg_gap.insertKmer(kmer_cauda);
                // buscando os kmers do gap
                //cout << "entre " << kmer_cabeca << " " << kmer_cauda << endl;
                for(int j = posicoesValidas[i] + 1; j < posicoesValidas[i+1]; j++)
                {
                    string kmer_aux = sequence.substr(j, k);
                    kmer_lista_aux = h.compareKmersWithGraph(kmer_aux, x * k); 
                    for (auto aux : kmer_lista_aux)
                        dbg_gap.insertKmer(aux);
                } 
                auto sequence_graph_aux = dbg_gap.dbgToSequenceGraph_1();        
                string sequence_aux = sequence.substr(posicoesValidas[i], k + dif);
                auto grafo_multicamada = m.buildMultilayerGraph(sequence_graph_aux.second, sequence_aux);
                auto retorno = m.dijkstra(grafo_multicamada, m.getInitialNode(), m.getEndNode());
                auto mapeamento = m.mostraMapeamento(retorno.first, sequence_graph_aux.first, sequence_graph_aux.second);    
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
            length = k * x, aux = 0;
            string kmer_cabeca = sequence.substr(posicoesValidas[posicoesValidas.size() - 1], k);
            dbg_gap.insertKmer(kmer_cabeca);   
            //kmers_mapeamento.push_back(kmer_cabeca);

            // buscando os kmers do gap
            for(int j = posicoesValidas[posicoesValidas.size() - 1]; j < sequence.length() - k; j++)
            {
                string kmer_aux = sequence.substr(j, k);
                kmer_lista_aux = h.compareKmersWithGraph(kmer_aux, x * k); 
                for (auto aux : kmer_lista_aux)
                    dbg_gap.insertKmer(aux);
            }  
            auto sequence_graph_aux = dbg_gap.dbgToSequenceGraph_1();
   
               int dif = sequence.length() - posicoesValidas[posicoesValidas.size() - 1];
            if (dif > 1)
            {    
                string sequence_aux = sequence.substr(posicoesValidas[posicoesValidas.size() - 1], k + dif);
                auto grafo_multicamada = m.buildMultilayerGraph(sequence_graph_aux.second, sequence_aux);
                auto retorno = m.dijkstra(grafo_multicamada, m.getInitialNode(), m.getEndNode());
                auto mapeamento = m.mostraMapeamento(retorno.first, sequence_graph_aux.first, sequence_graph_aux.second);    
                if (mapeamento.second.length() >= k)      
                    resposta = resposta + mapeamento.second.substr(k, mapeamento.second.length() - k);  
            }
        }
    } else
    {
        cout << "Nenhum kmer encontrado no grafo. A sequência vai ser toda mapeada." << endl;
        auto sequence_graph = h.dbgToSequenceGraph_1();
        auto grafo_multicamada = m.buildMultilayerGraph(sequence_graph.second, sequence);
        auto retorno = m.dijkstra(grafo_multicamada, m.getInitialNode(), m.getEndNode());
        auto resposta = m.mostraMapeamento(retorno.first, sequence_graph.first, sequence_graph.second);    
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
    // auto grafoSequencias = h.dbgToSequenceGraph_1();    
    // mapeamento
    auto retorno = mapeamento(h, sequence, k); 
    cout << "Mapeamento " << retorno << endl;
    return 0;
}