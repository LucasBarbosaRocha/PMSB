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
    int is_visited;
    string originalKmer;

    My_object();
};

My_object::My_object()
{
    this->adjacent = {};
    this->node_in_sequence_graph = 0;
    this->is_kmer_special = 0;
    this->node_kmer_for_kmer_special = 0;
    this->is_visited = 0;
}

class Hash
{
private:
    int k; /* comprimento do k-mer */
    unordered_map<string, My_object> new_graph; /* buckets */
    unordered_map<int, string> map_node_kmer;
    vector<string> oficialBases;
public:
    SequenceGraph sequenceGraph;
    SequenceGraph sequenceGraphReverse;
    /*  funcao construtor do grafo */
	Hash(int k); 

    string getKmerBySequenceGraphNode(int node);

    /* insere k-mers especiais no grafo */
    void enumerateKmers();

    vector<vector<double>> getCostMatrix(vector<int> ancoragem, string sequence, int k);

    /* insere k-mers especiais no grafo */
    void insertSpecialsKmers();

    void insertKmerInDbgAux(string kmer, My_object object);

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

    bool containsAndMarkIfExists(string kmer);

    /*  funcao transforma o grafo de De Bruijn em um grafo de sequências simples (um caractere por rótulo) */
    void dbgToTraditionalSequenceGraph(int reverse, int heuristic);

    /* funcao transforma o grafo de De Bruijn em um grafo de sequências simples reduzido */
    void dbgToSimplifiedSequenceGraph(int reverse);

    /* funcao imprime um grafo de De Bruijn */
	void displayHash();

    /* excluir grafo */
	void deleteHash();

    void deleteKmer(string kmer);

    string findKmerBySpecialNode(int index);

    string findKmerBySpecialKmer(string specialkmer);

    /* funcao devolve o comprimento k  */
    int getK();

    /* funcao recebe um arquivo com sequencias e insere todos os k-mers no grafo */
    void populateGraph(string nomeArquivo, bool detalhes);

    string readSequence(string nomeArquivo, bool detalhes);

    int getQtdKmers();

    int hammingDistance(string kmer1, string kmer2);

    void compareKmersWithGraph(Hash &dbg_gab, string kmer, int errors);

    string findKmerInTheGraph(string kmer, int errors);

    int compareKmersWithGraphAndRemove(Hash *dbg_aux, string kmer, int cost);

    void compareGraphWithSequence(Hash &dbg_gab, string sequence, int errors);

    string findKmerByIndex(int index);

    pair<int, int> vertice_inicial_final(string cabeca, string cauda);

    int verifyLabelExistisAndReturnVerticeIndex(string label);

    void insertKmersByNodes(list<int> nodes, Hash &dbg);

    int isKmerVisited(string kmer);

    int countVisitedKmers();

};

Hash::Hash(int k)
{
    this->k = k;
    this->new_graph = {};
    this->map_node_kmer = {};
    this->oficialBases = {"A", "C", "G", "T"};
}

string Hash::getKmerBySequenceGraphNode(int node)
{
    return this->map_node_kmer[node];
}

void Hash::enumerateKmers()
{
    int i = 0;
    for (auto it = this->new_graph.begin(); it != this->new_graph.end(); it++)
    {
        it->second.node_in_sequence_graph = i++;
    }
}

vector<vector<double>> Hash::getCostMatrix(vector<int> ancoragem, string sequence, int k)
{
    vector<vector<double>> matrix;
    int n = this->getQtdKmers();
    for (int line = 0; line < sequence.length() - (k-1); line++)
    {
        if (find(ancoragem.begin(), ancoragem.end(), line) == ancoragem.end())
        {   
            vector<double> v1;
            string kmer = sequence.substr(line, k);

            for (auto col = this->new_graph.begin(); col != this->new_graph.end(); col++)
            {
                cout << "buscando " << kmer << " " << line << " " << col->second.node_in_sequence_graph << " " << col->second.is_visited << " kmer " << col->first << endl;
                if (line != col->second.node_in_sequence_graph && col->second.is_visited != 1)
                {
                    int cost = this->hammingDistance(kmer, col->first);
                    v1.push_back(cost);
                } else
                    v1.push_back(3000);
            }
            matrix.push_back(v1);
            
            for (auto a : v1)
                cout << a << " ";
            cout << endl;
        }

    }
    return matrix;
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
    }    
}

void Hash::insertKmerInDbgAux(string kmer, My_object object)
{
    this->new_graph[kmer] = object;
}

void Hash::insertKmer(string kmer)
{
    list<string> adj; list<string>::iterator it; 
    My_object obj;
    string bases[] = {"A", "C", "G", "T"}, base_aux, kmer_aux;
    transform(kmer.begin(), kmer.end(),kmer.begin(), ::toupper);
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
                auto a = this->new_graph[kmer];
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
        if (kmer.find("N") == string::npos)
        {
            insertKmer(kmer);
        }
    }
}

bool Hash::containsAndMarkIfExists(string kmer)
{
    if (this->new_graph.find(kmer) != this->new_graph.end())
    {
        this->new_graph[kmer].is_visited = 1;
        return true;
    }
    //this->new_graph[kmer].is_visited = 0;
    return false;
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

void Hash::dbgToTraditionalSequenceGraph(int reverse, int heuristic)
{
    int qtdNodes = 0;
    unordered_map<string, My_object>:: iterator itr;
    string bases[] = {"A", "C", "G", "T"};
    unordered_map<int, string> kmerAndNode;

    for (auto itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        this->new_graph[itr->first].node_in_sequence_graph = qtdNodes;
        if (reverse == 0 && heuristic == 1)
            this->map_node_kmer[qtdNodes] = itr->first;
        qtdNodes += this->k;
    }

    if (reverse == 0)
        sequenceGraph.initilizeSequenceGraph(qtdNodes, this->k);
    else
        sequenceGraphReverse.initilizeSequenceGraph(qtdNodes, this->k);
    //SequenceGraph sequenceGraph(qtdNodes, this->k);

    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        int node = this->new_graph[itr->first].node_in_sequence_graph;
        for (int i = 0; i < itr->first.length(); i++)
        {
            if (i == 0)
            {
                if (reverse == 0)
                    sequenceGraph.alterarValorVerticeInicial(node + i, 1); // nao lembro para que fiz isso
                else
                    sequenceGraphReverse.alterarValorVerticeInicial(node + i, 1); // nao lembro para que fiz isso
            }
            if(reverse == 0)
                sequenceGraph.insertNode(node + i, itr->first.substr(i, 1)); 
            else
                sequenceGraphReverse.insertNode(node + i, itr->first.substr(i, 1)); 
          
            
            // kmerAndNode[aux.first + i] = aux.second; // mapeando kmer e node
        }
        for (int i = 0; i < itr->first.length() - 1; i++)
        {   if (reverse == 0)
                sequenceGraph.insertEdge(node + i, node + i + 1, 0);
            else
                sequenceGraphReverse.insertEdge(node + i + 1, node + i, 0);
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
                    sequenceGraphReverse.insertEdge(target_node, source_node, 0);
            }
        }   
    }
}


void Hash::dbgToSimplifiedSequenceGraph(int reverse)
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
    //SequenceGraph sequenceGraph(qtdNodes, this->k);
    if (reverse == 0)
        sequenceGraph.initilizeSequenceGraph(qtdNodes, this->k);
    else
        sequenceGraphReverse.initilizeSequenceGraph(qtdNodes, this->k);

    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        int node = this->new_graph[itr->first].node_in_sequence_graph;
        if (reverse == 0)
            sequenceGraph.insertNode(node, itr->first.substr(this->k-1,1));
        else
            sequenceGraphReverse.insertNode(node, itr->first.substr(this->k-1,1));
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
                    sequenceGraphReverse.insertEdge(target_node, source_node, 0);
            }
        }   
    }

    if (reverse == 0)
        sequenceGraph.markInitials(1);
    else
        sequenceGraphReverse.markInitials(1);
}

string Hash::findKmerBySpecialKmer(string specialkmer)
{
    return this->new_graph[specialkmer].originalKmer;
}

void Hash::deleteHash()
{
    this->new_graph.clear();
}

void Hash::deleteKmer(string kmer)
{
    this->new_graph.erase(kmer);
}

void Hash::populateGraph(string nomeArquivo, bool detalhes = false)
{
    string linha;
	fstream meuArquivo;
	meuArquivo.open(nomeArquivo);
    int qtd = 0;

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
        cout << "Qtd G.Kmers " << this->getQtdKmers() << endl;
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

void Hash::compareKmersWithGraph(Hash &dbg_gab, string kmer, int errors)
{
    unordered_map<string, My_object>:: iterator itr;
    list<string> kmers;

    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        if (itr->first.find("$") == string::npos)
        {
            if (hammingDistance(kmer, itr->first) < errors)
            {
                //cout << kmer << " " << itr->first << ":" << hammingDistance(kmer, itr->first) << " < " << errors << endl;
                dbg_gab.insertKmer(itr->first);
            }
        }
    }
}

string Hash::findKmerInTheGraph(string kmer, int errors)
{
    unordered_map<string, My_object>:: iterator itr;
    int cost = errors;
    list<string> kmers; string kmer_key = ""; int h = -1;

    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        if (itr->first.find("$") == string::npos)
        {
            h = hammingDistance(kmer, itr->first);
            if (h <= cost)
            {
                cost = h;
                kmer_key = itr->first;
            }
        }
    }
    return kmer_key;
}

int Hash::compareKmersWithGraphAndRemove(Hash *dbg_aux, string kmer, int cost)
{
    unordered_map<string, My_object>:: iterator itr;
    list<string> kmers; string kmer_key = ""; int h = -1;

    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        if (itr->first.find("$") == string::npos)
        {
            h = hammingDistance(kmer, itr->first);
            if (h <= cost)
            {
                cost = h;
                kmer_key = itr->first;
            }
        }
    }
    if (h != -1)
        this->deleteKmer(kmer_key);
    dbg_aux->insertKmer(kmer);
    return h;
}

void Hash::compareGraphWithSequence(Hash &dbg_gab, string sequence, int errors)
{
    unordered_map<string, My_object>:: iterator itr;
    list<string> kmers;
    string kmer; int val = 0;

    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        if (itr->first.find("$") == string::npos)
        {
            val = 0;
            for (int i = 0; i < sequence.size(); i+=k)
            {
                kmer = sequence.substr(i, k);          
                if (hammingDistance(kmer, itr->first) > errors)
                {
                    // cout << kmer << " " << itr->first << ":" << hammingDistance(kmer, itr->first) << " < " << errors << endl;
                    val = 1;
                }
            }
            if (val == 0)
            {
                dbg_gab.insertKmer(itr->first);
            }
        }
    }
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

string Hash::findKmerBySpecialNode(int index)
{
    unordered_map<string, My_object>:: iterator itr;
    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        if (itr->second.node_in_sequence_graph == index)
        {
            cout << "achei " << itr->second.node_in_sequence_graph  << " " << index << endl;
            return itr->first;
        }
    }
    return "";
}

pair<int, int> Hash::vertice_inicial_final(string cabeca, string cauda)
{
	int c = 0, t = 0; 
    unordered_map<string, My_object>:: iterator itr;

    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
	{
		if (itr->first.compare(cabeca) == 0)
        {
			c = itr->second.node_in_sequence_graph;
        }
		if (itr->first.compare(cauda) == 0)
        {
			t = itr->second.node_in_sequence_graph;
        }
	}
	return make_pair(c,t);
}

int Hash::verifyLabelExistisAndReturnVerticeIndex(string label)
{
    if (this->new_graph.find(label) != this->new_graph.end())
    {
        return this->new_graph[label].node_in_sequence_graph;
    }
    return -1;
}

void Hash::insertKmersByNodes(list<int> nodes, Hash &dbg)
{
    for (auto node : nodes)
    {
        for (int i = node - k; i < node + k; i++) 
        {
            if (i > 0) {
                string kmer = this->getKmerBySequenceGraphNode(i);
                if (kmer != "")
                {
                    if (this->new_graph[kmer].node_in_sequence_graph <= node && node <= this->new_graph[kmer].node_in_sequence_graph + this->k)
                    {
                        dbg.insertKmer(kmer);
                    }
                    break;
                } 
            }
        }
    } 
}

int Hash::isKmerVisited(string kmer)
{
    if (this->new_graph[kmer].is_visited == 1)
        return 1;
    return 0;
}

int Hash::countVisitedKmers()
{
    unordered_map<string, My_object>:: iterator itr;
    int count = 0;
    for (itr = this->new_graph.begin(); itr != this->new_graph.end(); itr++)
    {
        if(itr->second.is_visited == 0)
            count++;
    }
    return count;
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