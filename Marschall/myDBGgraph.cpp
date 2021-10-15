/*
    Classe myDBGgraph
    Autor: Lucas B. Rocha
    Orientador(es): Said Sadique e Francisco Elói
    Ano: 2021
    FACOM: Doutorado em Ciência da Computação    
    Implementação: 
        construção de um grafo de De Bruijnutilizando uma hash map 
*/
#include <bits/stdc++.h>
#include <math.h>       /* pow */
#include <string>
#include <unordered_map>
#include "mySequenceGraph.cpp"
#include <bits/stdc++.h>

using namespace std;

class Hash
{
private:
    int k; /* comprimento do k-mer */
    unordered_map<string, list<string>> graph; /* buckets */
    unordered_map<string, string> kmerSpecialAndKmer;
public:
    /*  funcao construtor do grafo */
	Hash(int k); 

    /* insere k-mers especiais no grafo */
    void insertSpecialsKmers();

    /* funcao recebe um k-mer e insere no grafo */
	void insertKmer(string kmer);

    /* funcao recebe uma sequencia e insere todos os k-mers no grafo  */
	void insertSequence(string sequence);    

    /* funcao recebe um k-mer kmer e uma base e verifica se tem aresta entre kmer e kmer[2,k] + base */
    bool containsOut(string kmer, string base); 

    /* funcao recebe um k-mer kamer e uma base e verifica se tem areste entre base + kmer[1,k-1] e kmer */
    bool containsIn(string kmer, string base);

    /* funcao recebe um k-mer kmer e verrifica se estah presente no grafo */
    bool contains(string kmer);

    /* funcao recebe um k-mer kmer e apaga do grafo  */
    void findAndRemoveKmer(string kmer);

    /*  funcao transforma o grafo de De Bruijn em um grafo de sequências simples (um caractere por rótulo) */
    pair<unordered_map<int, string>,SequenceGraph> dbgToSequenceGraph_1();

    /* funcao transforma o grafo de De Bruijn em um grafo de sequências simples reduzido */
    pair<unordered_map<int, string>,SequenceGraph> dbgToSequenceGraph_2();

    /* funcao imprime um grafo de De Bruijn */
	void displayHash();

    /* funcao devolve o comprimento k  */
    int getK();

    unordered_map<string, string> getKmerSpecialAndKmer();

    /* funcao recebe um arquivo com sequencias e insere todos os k-mers no grafo */
    void populateGraph(string nomeArquivo, bool detalhes);

    int getQtdKmers();

};

Hash::Hash(int k)
{
    this->k = k;
    graph = {};
}

void Hash::insertSpecialsKmers()
{
    unordered_map<string, list<string>>:: iterator itr;
    string bases[] = {"A", "C", "G", "T"}, kmer, kmer_aux; int val;
    list<string>::iterator it; list<string> bkp;

    for (itr = graph.begin(); itr != graph.end(); itr++)
        bkp.push_back((*itr).first);

    for (it = bkp.begin(); it != bkp.end(); it++)
    {
        kmer = *it; val = 1;
        for (auto base : bases)
        {
            //cout << "teste " << kmer.substr(k-1,1) << " " << base << " r " << kmer.substr(k-1,1).compare(base)<< endl;
            if (containsIn(kmer, base))
            {
                val = 0;
                break;
            }
        }
        if (val == 1)
        {
            //cout << " vamos inserrir " << kmer << endl;
            for (int i = 1; i < kmer.length(); i++)
            {
                kmer_aux = "";
                for (int j = 0; j < i; j++)
                    kmer_aux = kmer_aux + "$";                
                kmer_aux = kmer_aux + kmer.substr(0,this->k-i);
                kmerSpecialAndKmer[kmer_aux] = kmer;
                insertKmer(kmer_aux);
            }
        }
        else if(val == 1 && kmer.substr(0,this->k-1).compare(kmer.substr(1,this->k-1)) == 0)
        {
            for (int i = 1; i < kmer.length(); i++)
            {
                kmer_aux = "";
                for (int j = 0; j < i; j++)
                    kmer_aux = kmer_aux + "$";                
                kmer_aux = kmer_aux + kmer.substr(0,this->k-i);
                insertKmer(kmer_aux);
            }
        }

    }    
}

void Hash::insertKmer(string kmer)
{
    list<string> adj; list<string>::iterator it; 
    string bases[] = {"A", "C", "G", "T"}, base_aux, kmer_aux;
    if(!contains(kmer))
    {
        graph[kmer] = adj;
        for (auto base : bases)
        {
            if (containsOut(kmer, base))
            {             
                for (it = graph[kmer].begin(); it != graph[kmer].end(); it ++)
                    if (*it == base)
                        break;
                
                if (it == graph[kmer].end())
                    graph[kmer].push_back(base);
            }

            if (containsIn(kmer, base))
            {
                base_aux = kmer.substr(this->k-1,1); kmer_aux = base+kmer.substr(0,this->k-1);
                for (it = graph[kmer_aux].begin(); it != graph[kmer_aux].end(); it ++)
                    if (*it == base_aux)
                        break;
                
                if (it == graph[kmer_aux].end())
                    graph[kmer_aux].push_back(base_aux);
            }
        }
    } 
}

void Hash::insertSequence(string sequence)
{
    int n = sequence.length();
    string kmer;
    for (int i = 0; i <= n - k; i++)
    {
        kmer = sequence.substr(i, k);
        //cout << "kmer " << kmer << endl;
        insertKmer(kmer);
    }
}

bool Hash::contains(string kmer)
{
    if (graph.find(kmer) != graph.end())
        return true;
    return false;
}

bool Hash::containsOut(string kmer, string base)
{
    string kmer_aux = kmer.substr(1,(this->k-1)) + base;


    if (graph.find(kmer_aux) != graph.end())
        return true;
    return false;

}

bool Hash::containsIn(string kmer, string base)
{
    string kmer_aux = base + kmer.substr(0,(this->k-1));

    if (graph.find(kmer_aux) != graph.end())
        return true;
    return false;
}

pair<unordered_map<int, string>,SequenceGraph> Hash::dbgToSequenceGraph_1()
{
    int qtdNodes = 0;
    list<pair<int, string>> list_aux;
    unordered_map<string, list<string>>:: iterator itr;
    string bases[] = {"A", "C", "G", "T"};
    unordered_map<int, string> kmerAndNode;

    for (itr = graph.begin(); itr != graph.end(); itr++)
    {
        list_aux.push_back(make_pair(qtdNodes, itr->first));
        qtdNodes += this->k;
    }
    SequenceGraph sequenceGraph(qtdNodes, this->k);
    for (auto aux : list_aux)
    {
        for (int i = 0; i < aux.second.length(); i++)
        {
            if (i == 0)
            {
                sequenceGraph.alterarValorVerticeInicial(aux.first + i, 1);
            }
            sequenceGraph.insertNode(aux.first + i, aux.second.substr(i, 1)); 
            kmerAndNode[aux.first + i] = aux.second; // mapeando kmer e node
        }
        for (int i = 0; i < aux.second.length() - 1; i++)
            sequenceGraph.insertEdge(aux.first + i, aux.first + i + 1, 0);
    }

    for (auto aux : list_aux)
    {
        for (auto base : bases)
        {
            if(containsOut(aux.second, base))
            {
                int node = -1; string kmer_aux = aux.second.substr(1, this->k-1) + base;
                for (auto aux2 : list_aux)
                {
                    if (kmer_aux == aux2.second)
                    {
                        node = aux2.first + this -> k - 1;
                        break;
                    }
                }
                sequenceGraph.insertEdge(aux.first + this -> k - 1, node, 0);
            }
        }     
    }

    /*for (auto &a : kmerAndNode)
    {
        cout << a.first << " - " << a.second << endl;
    }*/

    return make_pair(kmerAndNode,sequenceGraph);
}


pair<unordered_map<int, string>,SequenceGraph> Hash::dbgToSequenceGraph_2()
{
    int qtdNodes = 0;
    list<pair<int, string>> list_aux;
    unordered_map<string, list<string>>:: iterator itr;
    string bases[] = {"A", "C", "G", "T"};
    unordered_map<int, string> kmerAndNode;

    this->insertSpecialsKmers();

    for (itr = graph.begin(); itr != graph.end(); itr++)
    {
        list_aux.push_back(make_pair(qtdNodes, itr->first));
        qtdNodes++;
    }
    SequenceGraph sequenceGraph(qtdNodes, this->k);

    for (auto aux : list_aux)
    {
        sequenceGraph.insertNode(aux.first, aux.second.substr(this->k-1,1));
        kmerAndNode[aux.first] = aux.second; // mapeando kmer e node
    }

    for (auto aux : list_aux)
    {
        for (auto base : bases)
        {
            if(containsOut(aux.second, base))
            {
                int node = -1; string kmer_aux = aux.second.substr(1, this->k-1) + base;
                for (auto aux2 : list_aux)
                {
                    if (kmer_aux == aux2.second)
                    {
                        node = aux2.first;
                        break;
                    }
                }
                sequenceGraph.insertEdge(aux.first, node, 0);
            }
        }     
    }

    sequenceGraph.markInitials();
    return make_pair(kmerAndNode,sequenceGraph);
}

void Hash::displayHash()
{
    unordered_map<string, list<string>>:: iterator itr;
    cout << "SequenceGraph De Bruijn: \n";
    list<string>::iterator it;
    for (itr = graph.begin(); itr != graph.end(); itr++)
    {
        cout << itr->first << ": ";
        for(it = itr->second.begin(); it != itr->second.end(); it++)
            cout << *it << " ";
        cout << endl;
    }
}

void Hash::populateGraph(string nomeArquivo, bool detalhes = false)
{
    string linha;
	fstream meuArquivo;
	meuArquivo.open(nomeArquivo);

	if (!meuArquivo) {
		cout << "Arquivo " << nomeArquivo << " de kmers não encontrado" << endl;
	}
	else {
        if (detalhes) cout << "Criando o grafo de De Bruijn" << endl;

        while (getline(meuArquivo, linha))
        {
            if (detalhes) cout << "Adicionando " << linha << endl;
            insertSequence(linha);
        }
		meuArquivo.close();
        if (detalhes) cout << "De Bruijn criado." << endl;
    }    
}


unordered_map<string, string> Hash::getKmerSpecialAndKmer()
{
    return this->kmerSpecialAndKmer;
}


int Hash::getQtdKmers()
{
    return this->graph.size();
}

/*int main()
{
    // array that contains kmers to be mapped
    string kmers[] = {"AAA", "AAT", "GGG", "GGT", "TTT", "ATT", "TTA", "TTT"};
    string kmers2[] = {"TACT","TTAA","ACGT","CGTG","GTTA","GTCT","TCAA", "TTT"};
    int n = sizeof(kmers2)/sizeof(kmers[0]);
    int k = 3;

    // insert the kmers into the hash table
    Hash h(k); // 7 is count of buckets in
                // hash table
    
    string nomeArquivo = "kmers.txt";

    h.populateGraph(nomeArquivo, false); 
    //h.insertSpecialsKmers();        
    h.displayHash();
    auto SequenceGraph = h.dbgToSequenceGraph_1();
    SequenceGraph.imprimeSequenceGraph();
return 0;
}*/