#include <iostream>
#include <queue>
#include <bits/stdc++.h>
#include <stack>
#include <list>
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"

vector<int> ancoras(Hash dbg, string sequence, int k)
{
    vector<int> ancoragem;
    for (int i = 0; i < sequence.length() - (k - 1); i++)
    {
        if (dbg.contains(sequence.substr(i,k)))
        {
            ancoragem.push_back(i);
        }
    }
    return ancoragem;
}

void mapeamento(Hash dbg, string sequence, int k, float x)
{
    int length, aux, dif = k, caminho_encontrado = 0;
    stack<string> pilha;
    list<string> kmers_mapeamento, kmer_lista_aux;
    string bases[] = {"A", "C", "G", "T"};
    Marschall m;
    Hash dbg_internal(k);
    int qtdKmers = dbg.getQtdKmers();

    for (int j = 0; j < sequence.size() - k; j++)
    {
        string kmer_aux = sequence.substr(j, k);
        kmer_lista_aux = dbg.compareKmersWithGraph(kmer_aux, x * k); 
        for (auto aux : kmer_lista_aux) {
            dbg_internal.insertKmer(aux);
        }
        kmer_lista_aux.clear();
    }

    auto sequence_graph_aux = dbg_internal.dbgToTraditionalSequenceGraph(0);    
    auto grafo_multicamada = m.buildMultilayerGraph(sequence_graph_aux, sequence);

    auto retorno = m.dijkstra(grafo_multicamada, m.getInitialNode(), m.getEndNode());
    auto mapeamento = m.showTraditionalMapping(retorno.first, dbg, sequence_graph_aux);

    if (retorno.second == INT_MAX)
        cout << "caminho nao encontrado" << endl;
    else
        cout << mapeamento.second << endl;
}

int main(int argc, char *argv[])
{
    MyUtils utils;

    if(utils.verifyData(argc, argv) == 1)
        exit (0);
    
    Hash h(utils.k);   
    utils.readSequence(utils.nameSequenceArchive);
    h.populateGraph(utils.nameArchive, false);
    float x = 0.5;
    mapeamento(h, utils.sequence, utils.k, x);
    return 0;
}