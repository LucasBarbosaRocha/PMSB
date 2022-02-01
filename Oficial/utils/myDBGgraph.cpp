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
#include <cmath>        // std::abs

using namespace std;

class My_object
{
public:
    list<string> adjacent;
    int node_in_sequence_graph;
    int is_kmer_special;
    int node_kmer_for_kmer_special;
    string originalKmer;

    My_object();
};

My_object::My_object()
{
    this->adjacent = {};
    this->node_in_sequence_graph = 0;
    this->is_kmer_special = 0;
    this->node_kmer_for_kmer_special = 0;
}

class Hash
{
private:
    int k; /* comprimento do k-mer */
    unordered_map<string, My_object> new_graph; /* buckets */
    vector<string> oficialBases;
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

    /*  funcao transforma o grafo de De Bruijn em um grafo de sequências simples (um caractere por rótulo) */
    SequenceGraph dbgToTraditionalSequenceGraph(int reverse);

    /* funcao transforma o grafo de De Bruijn em um grafo de sequências simples reduzido */
    SequenceGraph dbgToSimplifiedSequenceGraph(int reverse);

    /* funcao imprime um grafo de De Bruijn */
	void displayHash();

    /* excluir grafo */
	void deleteHash();

    string findKmerBySpecialKmer(string specialkmer);

    /* funcao devolve o comprimento k  */
    int getK();

    /* funcao recebe um arquivo com sequencias e insere todos os k-mers no grafo */
    void populateGraph(string nomeArquivo, bool detalhes);

    string readSequence(string nomeArquivo, bool detalhes);

    int getQtdKmers();

    int hammingDistance(string kmer1, string kmer2);

    list<string> compareKmersWithGraph(string kmer, int errors);

    string findKmerByIndex(int index);

    pair<int, int> vertice_inicial_final(string cabeca, string cauda);

};

Hash::Hash(int k)
{
    this->k = k;
    this->new_graph = {};
    this->oficialBases = {"A", "C", "G", "T"};
}

void Hash::insertSpecialsKmers()
{
    unordered_map<string, My_object>:: iterator itr;
    string bases[] = {"A", "C", "G", "T"}, kmer, kmer_aux; int val;

    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        kmer = (itr)->first; val = 1;
        for (auto base : bases)
        {
            //cout << "teste " << kmer.substr(k-1,1) << " " << base << " r " << kmer.substr(k-1,1).compare(base)<< endl;
            if (containsIn(kmer, base))
            {
                val = 0;
                break;
            }
        }
        if (val == 1 || (val == 1 && kmer.substr(0,this->k-1).compare(kmer.substr(1,this->k-1)) == 0))
        {
            this->new_graph[kmer].is_kmer_special = 1;
            for (int i = 1; i < kmer.length(); i++)
            {
                kmer_aux = "";
                for (int j = 0; j < i; j++)
                    kmer_aux = kmer_aux + "$";                
                kmer_aux = kmer_aux + kmer.substr(0,this->k-i);
                insertKmer(kmer_aux);
                this->new_graph[kmer_aux].is_kmer_special = 1;
                this->new_graph[kmer_aux].originalKmer = kmer;
            }
        }
        /*else if(val == 1 && kmer.substr(0,this->k-1).compare(kmer.substr(1,this->k-1)) == 0)
        {
            for (int i = 1; i < kmer.length(); i++)
            {
                kmer_aux = "";
                for (int j = 0; j < i; j++)
                    kmer_aux = kmer_aux + "$";                
                kmer_aux = kmer_aux + kmer.substr(0,this->k-i);
                insertKmer(kmer_aux);
            }
        }*/

    }    
}

void Hash::insertKmer(string kmer)
{
    list<string> adj; list<string>::iterator it; 
    My_object obj;
    string bases[] = {"A", "C", "G", "T"}, base_aux, kmer_aux;
    if(!contains(kmer))
    {
        obj.adjacent = adj; // new
        this->new_graph[kmer] = obj; // new
        for (auto base : bases)
        {
            if (containsOut(kmer, base))
            {             
                for (it = this->new_graph[kmer].adjacent.begin(); it != this->new_graph[kmer].adjacent.end(); it++)
                    if (*it == base)
                        break;
                
                if (it == this->new_graph[kmer].adjacent.end())
                    this->new_graph[kmer].adjacent.push_back(base);
            }

            if (containsIn(kmer, base))
            {
                base_aux = kmer.substr(this->k-1,1); kmer_aux = base+kmer.substr(0,this->k-1);
                for (it = this->new_graph[kmer].adjacent.begin(); it != this->new_graph[kmer].adjacent.end(); it++)
                    if (*it == base_aux)
                        break;
                
                if (it == this->new_graph[kmer].adjacent.end())
                    this->new_graph[kmer_aux].adjacent.push_back(base_aux);
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
    if (this->new_graph.find(kmer) != this->new_graph.end())
        return true;
    return false;
}

bool Hash::containsOut(string kmer, string base)
{
    string kmer_aux = kmer.substr(1,(this->k-1)) + base;
    if (this->new_graph.find(kmer_aux) != this->new_graph.end())
    {
        return true;
    }
    return false;

}

bool Hash::containsIn(string kmer, string base)
{
    string kmer_aux = base + kmer.substr(0,(this->k-1));

    if (this->new_graph.find(kmer_aux) != this->new_graph.end())
        return true;
    return false;
}

SequenceGraph Hash::dbgToTraditionalSequenceGraph(int reverse)
{
    int qtdNodes = 0;
    unordered_map<string, My_object>:: iterator itr;
    string bases[] = {"A", "C", "G", "T"};
    unordered_map<int, string> kmerAndNode;

    for (auto itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        this->new_graph[itr->first].node_in_sequence_graph = qtdNodes;
        qtdNodes += this->k;
    }
    SequenceGraph sequenceGraph(qtdNodes, this->k);

    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        int node = this->new_graph[itr->first].node_in_sequence_graph;
        for (int i = 0; i < itr->first.length(); i++)
        {
            if (i == 0)
            {
                sequenceGraph.alterarValorVerticeInicial(node + i, 1); // nao lembro para que fiz isso
            }
            sequenceGraph.insertNode(node + i, itr->first.substr(i, 1)); 
            // kmerAndNode[aux.first + i] = aux.second; // mapeando kmer e node
        }
        for (int i = 0; i < itr->first.length() - 1; i++)
        {   if (reverse == 0)
                sequenceGraph.insertEdge(node + i, node + i + 1, 0);
            else
                sequenceGraph.insertEdge(node + i + 1, node + i, 0);
        }
    }

    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        for (auto base : bases)
        {
            if(containsOut(itr->first, base))
            {
                string kmer_aux = itr->first.substr(1, this->k-1) + base;
                int source_node = this->new_graph[itr->first].node_in_sequence_graph + this -> k - 1;
                int target_node = this->new_graph[kmer_aux].node_in_sequence_graph + this -> k - 1;
                if (reverse == 0)
                    sequenceGraph.insertEdge(source_node, target_node, 0);
                else 
                    sequenceGraph.insertEdge(target_node, source_node, 0);
            }
        }   
    }
    return sequenceGraph;
}


SequenceGraph Hash::dbgToSimplifiedSequenceGraph(int reverse)
{
    int qtdNodes = 0;
    unordered_map<string, My_object>:: iterator itr;
    string bases[] = {"A", "C", "G", "T"};

    this->insertSpecialsKmers();

    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        this->new_graph[itr->first].node_in_sequence_graph = qtdNodes;
        qtdNodes++;
    }
    SequenceGraph sequenceGraph(qtdNodes, this->k);

    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        int node = this->new_graph[itr->first].node_in_sequence_graph;
        sequenceGraph.insertNode(node, itr->first.substr(this->k-1,1));
    }

    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        for (auto base : bases)
        {
            if(containsOut(itr->first, base))
            {
                string kmer_aux = itr->first.substr(1, this->k-1) + base;
                int source_node = this->new_graph[itr->first].node_in_sequence_graph;
                int target_node = this->new_graph[kmer_aux].node_in_sequence_graph;
                if (reverse == 0)
                    sequenceGraph.insertEdge(source_node, target_node, 0);
                else 
                    sequenceGraph.insertEdge(target_node, source_node, 0);
            }
        }   
    }
    sequenceGraph.markInitials(1);
    return sequenceGraph;
}

string Hash::findKmerBySpecialKmer(string specialkmer)
{
    return this->new_graph[specialkmer].originalKmer;
}

void Hash::deleteHash()
{
    this->new_graph.clear();
}

void Hash::populateGraph(string nomeArquivo, bool detalhes = false)
{
    string linha;
	fstream meuArquivo;
	meuArquivo.open(nomeArquivo);
    int qtd = 0;
    cout << "nome arquivo " << nomeArquivo << endl;
	if (!meuArquivo) {
		cout << "Arquivo " << nomeArquivo << " de kmers não encontrado" << endl;
	}
	else {
        if (detalhes) cout << "Criando o grafo de De Bruijn" << endl;
        getline(meuArquivo, linha);
        while (getline(meuArquivo, linha))
        {
            if (detalhes) cout << "Adicionando " << linha << endl;
            insertSequence(linha);
            qtd = qtd + linha.size() - k;
            getline(meuArquivo, linha);
        }
		meuArquivo.close();
        if (detalhes) cout << "De Bruijn criado." << endl;
        cout << "Qtd. Kmers " << this->getQtdKmers() << endl;
    }    
}

string Hash::readSequence(string nomeArquivo, bool detalhes = false)
{
    string linha;
	fstream meuArquivo;
	meuArquivo.open(nomeArquivo);

	if (!meuArquivo) {
		cout << "Arquivo " << nomeArquivo << " de kmers não encontrado" << endl;
	}
	else {
        if (detalhes) cout << "Criando o grafo de De Bruijn" << endl;

        getline(meuArquivo, linha);     
        getline(meuArquivo, linha);     
		meuArquivo.close();
        if (detalhes) cout << "De Bruijn criado." << endl;
    }    
    return linha;
}

int Hash::getQtdKmers()
{
    return this->new_graph.size();
}

int Hash::hammingDistance(string kmer1, string kmer2)
{
    int i, errors = 0;
    for (i = 0; i < kmer1.length(); i++)
    {
        if (kmer1[i] != kmer2[i])
        {
            errors++;
        }
    }
    return errors;
}

list<string> Hash::compareKmersWithGraph(string kmer, int errors)
{
    unordered_map<string, My_object>:: iterator itr;
    list<string> kmers;
    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        //cout << "Bla " << kmer << " " << itr->first << endl;
        if (hammingDistance(kmer, itr->first) <= errors)
        {
            kmers.push_front(itr->first);
        }
    }
    return kmers;
}

void Hash::displayHash()
{
    unordered_map<string, My_object>:: iterator itr;
    cout << "SequenceGraph De Bruijn: \n";
    list<string>::iterator it;
    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        cout << itr->first << ": ";
        for(it = itr->second.adjacent.begin(); it != itr->second.adjacent.end(); it++)
            cout << *it << " ";
        cout << endl;
    }
}

string Hash::findKmerByIndex(int index)
{
    unordered_map<string, My_object>:: iterator itr;
    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        int dif = abs(itr->second.node_in_sequence_graph - index);
        //cout << "no teste " << itr->first << " " << itr->second.node_in_sequence_graph  << " <=> " << index << "dif " << dif << endl;
        if (dif <= this->k)
        {
            //cout << "kmer teste " << itr->first << endl;
            return itr->first;
        }
    }
    return "";
}

pair<int, int> Hash::vertice_inicial_final(string cabeca, string cauda)
{
	list<int> head, tail; 
    unordered_map<string, My_object>:: iterator itr;

    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
	{
		if (itr->first.compare(cabeca))
			head.push_back(itr->second.node_in_sequence_graph);
		if (itr->first.compare(cauda))
			head.push_back(itr->second.node_in_sequence_graph);
	}
	head.sort();
	tail.sort();
	int c = head.front();
	tail.reverse();
	int t = tail.front();
	head.clear();
	tail.clear();
	return make_pair(c,t);
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
    
    string nomeArquivo = "../archives/kmers.txt";

    h.populateGraph(nomeArquivo, false); 
    //h.insertSpecialsKmers();        
    h.displayHash();
    auto SequenceGraph = h.dbgToTraditionalSequenceGraph();
    SequenceGraph.printGraph();

    auto SequenceGraph2 = h.dbgToSimplifiedSequenceGraph();
    h.displayHash();
    SequenceGraph2.printGraph(); 
return 0;
} */