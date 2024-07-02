// de Bruijn sequence mapping Tool with graph Node change

#include <fstream>
#include <iostream>
#include "utils/marschall.cpp"
#include "utils/myUtils.cpp"
#include <list>
using namespace std;

MyUtils utils;

pair<int, list<string>> exato(Hash h, string kmer_sequence, int k, bool detalhes)
{
    int qtd = 0;
    list<string> kmers;
    if (detalhes) cout << "Procurando kmers da sequência " << kmer_sequence << " no grafo" << endl;
    
    auto rep = h.findAnchors(kmer_sequence);
    cout << "Ancoras " << rep.size() << endl;

    for (int i = 0; i < kmer_sequence.length()-(k-1); i++)
    {
        string kmer = kmer_sequence.substr(i,k);
      
        if(!h.contains(kmer))
        {
            if (detalhes) cout << "Vamos add.: " << kmer << endl;
            kmers.push_back(kmer);
            h.insertSequence(kmer_sequence.substr(i,k));
            qtd++;
        }else
        {
            if (detalhes) cout << "Tem: " << kmer << endl;
        } 
    }

    rep = h.findAnchors(kmer_sequence);
    cout << "Ancoras " << rep.size() << endl;

    return make_pair(qtd, kmers);
}

int main(int argc, char *argv[])
{
    string line;
    if(utils.verifyData(argc, argv) == 1)
        exit (0); 
        
    Hash h(utils.k);   
    ifstream file(utils.nameSequenceArchive);
    h.populateGraph(utils.nameArchive, false);    

    while(getline(file, line))
    {
        //utils.readSequence(utils.nameSequenceArchive);  
        getline(file, line);
        cout << "Kmers: " << h.getQtdKmers() << endl;
        cout << "Size L.Read " << line.size() << endl;
        utils.sequence = line;
        // mapeamento
        auto retorno = exato(h, utils.sequence, utils.k, false);
        cout << "Precisamos inserir: " << retorno.first  << " kmers" << endl;
    }     
    return 0;
}



