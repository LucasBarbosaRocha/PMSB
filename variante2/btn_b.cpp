// de Bruijn sequence mapping Tool with graph Node change

#include <fstream>
#include <iostream>
#include "utils/myUtils.cpp"
#include "utils/myBifrost/my_bifrost.cpp"
#include <list>
using namespace std;

MyUtils utils;

pair<int, list<string>> exato(my_Bifrost bf, string kmer_sequence, int k, bool detalhes)
{
    int qtd = 0;
    list<string> kmers;
    if (detalhes) cout << "Procurando kmers da sequência " << kmer_sequence << " no grafo" << endl;

    auto rep = bf.findAnchors(kmer_sequence);
    cout << "Ancoras " << rep.size() << endl;

    if (kmer_sequence.length() < (size_t)k)
        return make_pair(qtd, kmers);
    for (size_t i = 0; i <= kmer_sequence.length() - k; i++)
    {
        const Kmer kmer = Kmer(kmer_sequence.substr(i,k).c_str());
      
        if(!bf.haskmer(kmer))
        {
            if (detalhes) cout << "Vamos add.: " << kmer.toString() << endl;
            kmers.push_back(kmer.toString());
            bf.insertSequence(kmer_sequence.substr(i,k));
            qtd++;
        }else
        {
            if (detalhes) cout << "Tem: " << kmer.toString() << endl;
        } 
    }

    rep = bf.findAnchors(kmer_sequence);
    cout << "Ancoras " << rep.size() << endl;

    return make_pair(qtd, kmers);
}

int main(int argc, char *argv[])
{
    string line;
    if(utils.verifyData(argc, argv) == 1)
        exit (0); 
        
    ifstream file(utils.nameSequenceArchive);
    my_Bifrost bf(utils.k, utils.nameArchive);  

    while(getline(file, line))
    {
        //utils.readSequence(utils.nameSequenceArchive);  
        getline(file, line);
        cout << "Kmers: " << bf.size() << endl;
        cout << "Size L.Read " << line.size() << endl;
        utils.sequence = line;
        // mapeamento
        auto retorno = exato(bf, utils.sequence, utils.k, false);
        cout << "Precisamos inserir: " << retorno.first  << " kmers" << endl;
    }     
    return 0;
}



