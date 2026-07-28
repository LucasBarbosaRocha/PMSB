#include <iostream>
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"
#include <queue>
#include <bits/stdc++.h>

using namespace std;

pair<int, list<string>> exato(Hash &cdbg, const string &kmer_sequence, int k, bool detalhes);

int main(int argc, char **argv)
{
    MyUtils utils;
    if (utils.verifyData(argc, argv) == 1)
        exit(0);

    Hash h(utils.k);
    utils.readSequence(utils.nameSequenceArchive);
    h.populateGraph(utils.nameArchive, false);

    if(utils.typeGraph == 0)
    {
        auto retorno = exato(h, utils.sequence, utils.k, false);
        cout << "Exato: " << endl;
        cout << "Vamos inserir: " << retorno.first << " ";
    }
    return 0;
}

pair<int, list<string>> exato(Hash &cdbg, const string &kmer_sequence, int k, bool detalhes)
{
    int qtd = 0;
    list<string> kmers;
    if (detalhes) cout << "Procurando kmers da sequência " << kmer_sequence << " no grafo" << endl;
    for (int i = 0; i < kmer_sequence.length()-(k-1); i++)
    {
        string kmer = kmer_sequence.substr(i,k);

        if(!cdbg.contains(kmer))
        {
            if (detalhes) cout << "Vamos add.: " << kmer << endl;
            kmers.push_back(kmer);
            qtd++;
        }else
        {
            if (detalhes) cout << "Tem: " << kmer << endl;
        } 
    }
    return make_pair(qtd, kmers);
}
