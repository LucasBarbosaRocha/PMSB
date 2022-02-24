#include <iostream>
#include <queue>
#include <bits/stdc++.h>
#include <stack>
#include <list>
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"

#define x 0.5
MyUtils utils;

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

string mapeamento(Hash dbg, string sequence, int k)
{
    int length, aux, dif = k, caminho_encontrado = 0;
    stack<string> pilha;
    list<string> kmers_mapeamento, kmer_lista_aux;
    string bases[] = {"A", "C", "G", "T"};
    Marschall m;
    Hash dbg_internal(k);
    pair<list<string>, string> mapeamento;
    int qtdKmers = dbg.getQtdKmers();

    for (int j = 0; j < sequence.size() - k; j++)
    {
        string kmer_aux = sequence.substr(j, k);
        dbg.compareKmersWithGraph(dbg_internal, kmer_aux, x * k); 
        for (auto aux : kmer_lista_aux) {
            dbg_internal.insertKmer(aux);
        }
        kmer_lista_aux.clear();
    }

    if (utils.typeGraph == 0)
        dbg_internal.dbgToTraditionalSequenceGraph(0);   
    else
        dbg_internal.dbgToSimplifiedSequenceGraph(0);    
    

    m.buildMultilayerGraph(dbg_internal.sequenceGraph, sequence);

    auto retorno = m.dijkstra(m.m_sequenceGraph, m.getInitialNode(), m.getEndNode());
    if (utils.typeGraph == 0)
        mapeamento = m.showTraditionalMapping(retorno.first, dbg, dbg_internal.sequenceGraph);
    else
        mapeamento = m.showTraditionalMapping(retorno.first, dbg, dbg_internal.sequenceGraph);

    if (mapeamento.second.size() > 0)
        return mapeamento.second;
    else
        return "Sequência nao pode ser mapeada";
}

int main(int argc, char *argv[])
{
    if(utils.verifyData(argc, argv) == 1)
        exit (0);
    
    Hash h(utils.k);   
    utils.readSequence(utils.nameSequenceArchive);
    h.populateGraph(utils.nameArchive, false);
    auto retorno = mapeamento(h, utils.sequence, utils.k);
    cout << retorno << endl;
    return 0;
}