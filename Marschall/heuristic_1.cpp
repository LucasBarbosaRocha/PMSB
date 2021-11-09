#include <iostream>
#include "myDBGgraph.cpp"
#include <queue>
#include <bits/stdc++.h>
#include <stack>
#include <list>

vector<int> ancoras(Hash dbg, string sequence, int k)
{
    vector<int> ancoragem;
    for (int i = 0; i < sequence.length() - (k - 1); i++)
    {
        if (dbg.contains(sequence.substr(i,k)))
        {
            ancoragem.push_back(i);
            //cout << "Ancora " << i << ": " << sequence.substr(i,k) << endl;
        }
    }
    return ancoragem;
}

void mapeamento(Hash dbg, string sequence, int k, float x)
{
    int length, aux, dif, caminho_encontrado = 0;
    stack<string> pilha;
    list<string> kmers_mapeamento, kmer_lista_aux;
    string bases[] = {"A", "C", "G", "T"};
    Hash dbg_gap(k);


    auto ancoragem = ancoras(dbg, sequence, k);
    auto sequence_graph_aux_3 = dbg_gap.dbgToSequenceGraph_3();        
    auto sequence_graph_aux_3R = dbg_gap.dbgToSequenceGraph_3R();        

    list<string> kmer_lista_oficial;
    // cabeca sem ancora

    if(ancoragem.size() > 0 && ancoragem[0] != 0)
    {
        length = dif * x, aux = 0;

        string kmer_cauda = sequence.substr(ancoragem[0], k);

        dbg_gap.insertKmer(kmer_cauda);   
        // buscando os kmers do gap

        for (int j = 0; j < sequence_graph_aux_3.getV(); j++)
        {
            if (sequence_graph_aux_3.isExistsPathLength(0, length) && sequence_graph_aux_3R.isExistsPathLength(0, length))
            {
                dbg_gap.insertKmer(sequence_graph_aux_3.getBase(j));
            }
        }

        auto sequence_graph_aux = dbg_gap.dbgToSequenceGraph_3();    
    
        for (auto it = kmer_lista_aux.begin(); it != kmer_lista_aux.end();)
        {
            int inicial = sequence_graph_aux.verifyLabelExistisAndReturnVerticeIndex(*it);
            int final = sequence_graph_aux.verifyLabelExistisAndReturnVerticeIndex(kmer_cauda);
            auto r = sequence_graph_aux.dijkstra(inicial, final);     
            if (r.second == INT_MAX)
            {
                it++;
            } else {
                int val = 0;
                for (auto r2 : r.first)
                {
                    if (val != 0)
                        kmers_mapeamento.push_front(r2.second);
                    val++;
                }
                caminho_encontrado = 1;
                it = kmer_lista_aux.end();
            }    
        }   
    }

    // ancoras
    if (ancoragem.size() > 0)
    {
        for(int i = 0; i < ancoragem.size() - 1; i++)
        {     

            dif = ancoragem[i + 1] -  ancoragem[i], caminho_encontrado = 0;
            kmer_lista_aux.clear();
            string kmer_cabeca = sequence.substr(ancoragem[i], k);
            kmers_mapeamento.push_back(kmer_cabeca);
            if (dif > 1)
            {
                length = dif * x, aux = 0;
                string kmer_cabeca = sequence.substr(ancoragem[i], k);
                pilha.push(kmer_cabeca);
                string kmer_cauda = sequence.substr(ancoragem[i+1], k);
                dbg_gap.deleteHash();
                dbg_gap.insertKmer(kmer_cabeca);
                dbg_gap.insertKmer(kmer_cauda);
                // buscando os kmers do gap
                //cout << "entre " << kmer_cabeca << " " << kmer_cauda << endl;
                for (int j = 0; j < sequence_graph_aux_3.getV(); j++)
                {
                    if (sequence_graph_aux_3.isExistsPathLength(0, length) && sequence_graph_aux_3R.isExistsPathLength(0, length))
                    {
                        dbg_gap.insertKmer(sequence_graph_aux_3.getBase(j));
                    }
                }

                auto sequence_graph_aux = dbg_gap.dbgToSequenceGraph_3();        
                int inicial = sequence_graph_aux.verifyLabelExistisAndReturnVerticeIndex(kmer_cabeca);
                int final = sequence_graph_aux.verifyLabelExistisAndReturnVerticeIndex(kmer_cauda);
                auto r = sequence_graph_aux.dijkstra(inicial, final);
                //sequence_graph_aux.printGraph();

                if (r.second == INT_MAX)
                {
                    //cout << "caminho nao encontrado " << endl;
                    caminho_encontrado = 1;
                    i = ancoragem.size() + 1;
                } else {
                    //cout << "caminho encontrado de custo " << r.second << endl;
                    int val = 0;
                    for (auto r2 : r.first)
                    {
                        if (val != 0)
                            kmers_mapeamento.push_front(r2.second);
                        val++;
                    }
                }
            }
        }  
    } 

    //cauda sem ancora no final
    if ( ancoragem.size() > 0 && caminho_encontrado != 1 && ancoragem[ancoragem.size() - 1] != sequence.length() - k)
    {
        dif = sequence.length() - ancoragem[ancoragem.size() - 1];
        dbg_gap.deleteHash();
        if (dif > 1)
        {
            length = dif * x, aux = 0;
            string kmer_cabeca = sequence.substr(ancoragem[ancoragem.size() - 1], k);
            dbg_gap.insertKmer(kmer_cabeca);   
            //kmers_mapeamento.push_back(kmer_cabeca);

            // buscando os kmers do gap
            for (int j = 0; j < sequence_graph_aux_3.getV(); j++)
            {
                if (sequence_graph_aux_3.isExistsPathLength(0, length) && sequence_graph_aux_3R.isExistsPathLength(0, length))
                {
                    dbg_gap.insertKmer(sequence_graph_aux_3.getBase(j));
                }
            } 

            auto sequence_graph_aux = dbg_gap.dbgToSequenceGraph_3();
            for (auto it = kmer_lista_aux.begin(); it != kmer_lista_aux.end();)
            {
                int inicial = sequence_graph_aux.verifyLabelExistisAndReturnVerticeIndex(kmer_cabeca);
                int final = sequence_graph_aux.verifyLabelExistisAndReturnVerticeIndex(*it);
                auto r = sequence_graph_aux.dijkstra(inicial, final);     
    
                if (r.second == INT_MAX)
                {
                    //cout << "nao achei caminho nao encontrado " << endl;
                    it++;
                } else {
                    //cout << "achei caminho encontrado de custo " << r.second << endl;
                    int val = 0;
                    for (auto r2 : r.first)
                    {
                        if (val != 0)
                            kmers_mapeamento.push_back(r2.second);
                        val++;
                    }
                    it = kmer_lista_aux.end();
                } 
            } 
        }      
    }

    if (kmers_mapeamento.size() == 0)
    {
        cout << "sequencia não pode ser mapeada :-/" << endl;
    }else
    {
        cout << "Mapeamento " << endl;
        int val = 0;
        for (auto a : kmers_mapeamento)
        {
            if (val == 0)
                cout << a;
            else 
                cout << a[a.length() - 1];
            val++;
        }
        cout << endl;
    }
}

int main (void)
{
    string sequence = "TGCAAAACTATAAAAATTTTGGAACGAGTCAAAGGGTAGGTCTTGGGAATACTACTTCATCCGGTTCTTG", nomeArquivo = "/home/lucas/Documentos/sequencias/sequencias100.fasta";
    // AAGACGTGTAGACGTGA
    map<string, int> kmersInTheGraph;
    int k = 10;
    float x = 0.5;
    Hash h(k);   
    h.populateGraph(nomeArquivo, false);
    auto grafoSequencias = h.dbgToSequenceGraph_3();    
    // mapeamento
    mapeamento(h, sequence, k, x);
    return 0;
}