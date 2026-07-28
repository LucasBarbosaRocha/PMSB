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

constexpr bool kVerboseDbg = false;

class My_object
{
public:
    list<string> adjacent;
    int node_in_sequence_graph;
    int is_kmer_special;
    int node_kmer_for_kmer_special;
    int is_visited;
    string originalKmer;
    int ordem;

    My_object();
};

My_object::My_object()
{
    this->adjacent = {};
    this->node_in_sequence_graph = 0;
    this->is_kmer_special = 0;
    this->node_kmer_for_kmer_special = 0;
    this->is_visited = 0;
    this->ordem = 0;
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
	void insertKmer(string kmer, int ordem);

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
                //cout << "buscando " << kmer << " " << line << " " << col->second.node_in_sequence_graph << " " << col->second.is_visited << " kmer " << col->first << endl;
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
    string bases[] = {"A", "C", "G", "T"};

    for (auto& pair : new_graph)
    {
        const string& kmer = pair.first;
        bool isSpecial = true;
        
        for (const auto& base : bases)
        {
            if (containsIn(kmer, base))
            {
                isSpecial = false;
                break;
            }
        }
        
        if (isSpecial || (isSpecial && kmer.substr(0, this->k - 1) == kmer.substr(1, this->k - 1)))
        {
            auto& node = new_graph[kmer];
            node.is_kmer_special = 1;
            int ordem = node.ordem;
            
            for (int i = 1; i < kmer.length(); ++i)
            {
                string kmer_aux(i, '$');
                kmer_aux += kmer.substr(0, this->k - i);
                int new_ordem = ordem - i;
                insertKmer(kmer_aux, new_ordem);
                auto& new_node = new_graph[kmer_aux];
                new_node.is_kmer_special = 1;
                new_node.originalKmer = kmer;
            }
        }
    }
}

void Hash::insertKmerInDbgAux(string kmer, My_object object)
{
    this->new_graph[kmer] = object;
}

void Hash::insertKmer(string kmer, int ordem) {
    string bases[] = {"A", "C", "G", "T"};
    string base_aux, kmer_aux;
    string upper_kmer = kmer;
    transform(upper_kmer.begin(), upper_kmer.end(), upper_kmer.begin(), ::toupper);

    if (!contains(upper_kmer)) {
        My_object obj;
        obj.ordem = ordem;
        new_graph[upper_kmer] = obj;

        auto& adjacent = new_graph[upper_kmer].adjacent;
        unordered_set<string> adj_set(adjacent.begin(), adjacent.end());  // Use a set to avoid duplicates

        for (const auto& base : bases) {
            if (containsOut(upper_kmer, base)) {
                if (adj_set.find(base) == adj_set.end()) {
                    adjacent.push_back(base);
                    adj_set.insert(base);
                }
            }

            if (containsIn(upper_kmer, base)) {
                base_aux = upper_kmer.substr(this->k - 1, 1);
                kmer_aux = base + upper_kmer.substr(0, this->k - 1);
                auto& kmer_aux_adjacent = new_graph[kmer_aux].adjacent;
                unordered_set<string> kmer_aux_adj_set(kmer_aux_adjacent.begin(), kmer_aux_adjacent.end());

                if (kmer_aux_adj_set.find(base_aux) == kmer_aux_adj_set.end()) {
                    kmer_aux_adjacent.push_back(base_aux);
                    kmer_aux_adj_set.insert(base_aux);
                }
            }
        }
    }
}

void Hash::insertSequence(string sequence)
{
    int n = sequence.length();
    string kmer;

    if (sequence.find("N") != string::npos) {
        return;
    }

    for (int i = 0; i <= n - k; ++i)
    {
        kmer = sequence.substr(i, k);
        insertKmer(kmer, i);
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
    string bases[] = {"A", "C", "G", "T"};
    unordered_map<int, string> kmerAndNode;
    auto& graph = (reverse == 0) ? sequenceGraph : sequenceGraphReverse;

    // Map nodes and optionally map k-mers
    for (auto& pair : new_graph) {
        auto& key = pair.first;
        pair.second.node_in_sequence_graph = qtdNodes;
        if (reverse == 0 && heuristic == 1) {
            map_node_kmer[qtdNodes] = key;
        }
        qtdNodes += this->k;
    }

    graph.initilizeSequenceGraph(qtdNodes, this->k);

    // Insert nodes and initial edges
    for (auto& pair : new_graph) {
        auto& key = pair.first;
        int node = pair.second.node_in_sequence_graph;
        for (int i = 0; i < this->k; ++i) {
            if (i == 0) {
                graph.alterarValorVerticeInicial(node + i, 1);
            }
            graph.insertNode(node + i, key.substr(i, 1));
            if (i < this->k - 1) {
                if (reverse == 0) {
                    graph.insertEdge(node + i, node + i + 1, 0);
                } else {
                    graph.insertEdge(node + i + 1, node + i, 0);
                }
            }
        }
    }

    // Insert remaining edges based on k-mer connections
    for (auto& pair : new_graph) {
        auto& key = pair.first;
        int source_node = pair.second.node_in_sequence_graph + this->k - 1;
        for (const auto& base : bases) {
            if (containsOut(key, base)) {
                string kmer_aux = key.substr(1, this->k - 1) + base;
                int target_node = new_graph[kmer_aux].node_in_sequence_graph + this->k - 1;
                if (reverse == 0) {
                    graph.insertEdge(source_node, target_node, 0);
                } else {
                    graph.insertEdge(target_node, source_node, 0);
                }
            }
        }
    }
}

void Hash::dbgToSimplifiedSequenceGraph(int reverse) {
    int qtdNodes = 0, novosNodes = 0;
    unordered_map<string, My_object>::iterator itr;
    string bases[] = {"A", "C", "G", "T"};

    this->insertSpecialsKmers();

    vector<pair<string, int>> keys;
    for (const auto& pair : new_graph) {
        keys.emplace_back(pair.first, pair.second.ordem);
    }

    sort(keys.begin(), keys.end(), [](const pair<string, int>& a, const pair<string, int>& b) {
        return a.second < b.second;
    });

    for (const auto& key_pair : keys) {
        const string& key = key_pair.first;
        auto& node = new_graph[key];
        //cout << key << "" << qtdNodes << endl;
        node.node_in_sequence_graph = qtdNodes;

        if (all_of(key.begin(), key.end(), [&](char c) { return c == key[0]; })) {
            qtdNodes++;
        }
        qtdNodes++;
    }
    //qtdNodes += novosNodes;

    auto& graph = (reverse == 0) ? sequenceGraph : sequenceGraphReverse;
    graph.initilizeSequenceGraph(qtdNodes, this->k);

    int nosInseridos = 0;
    for (const auto& pair : new_graph) {
        const string& key = pair.first;
        int node = pair.second.node_in_sequence_graph;
        graph.insertNode(node, key.substr(this->k-1, 1));
        nosInseridos++;
    }

    for (const auto& pair : new_graph) {
        const string& key = pair.first;
        int source_node = pair.second.node_in_sequence_graph;
        for (const auto& base : bases) {
            if (containsOut(key, base)) {
                string kmer_aux = key.substr(1, this->k-1) + base;
                int target_node = new_graph[kmer_aux].node_in_sequence_graph;

                if (reverse == 0) {
                    if (kVerboseDbg) cout << source_node << " -> " << target_node << endl;
                    if (source_node != target_node) {
                        graph.insertEdge(source_node, target_node, 0);
                        graph.insertIncoming(source_node, target_node);
                        graph.insertOutComing(source_node, target_node);
                    } else {
                        graph.insertNode(source_node+1, key.substr(this->k-1, 1));
                        graph.insertEdge(source_node, source_node+1, 0);
                        graph.insertEdge(source_node+1, source_node, 0);
                        graph.insertIncoming(source_node + 1, source_node);                  
                        graph.insertOutComing(source_node + 1, source_node);
                        nosInseridos++;
                    }
                } else {
                    graph.insertEdge(target_node, source_node, 0);
                }
            }
        }
    }

    for (const auto& pair : new_graph) {
        const string& key = pair.first;
        int source_node = pair.second.node_in_sequence_graph;
        for (const auto& base : bases) {
            if (containsOut(key, base)) {
                string kmer_aux = key.substr(1, this->k-1) + base;
                int target_node = new_graph[kmer_aux].node_in_sequence_graph;

                if (reverse == 0) {
                    if (kVerboseDbg) cout << source_node << " -> " << target_node << endl;
                    if (source_node == target_node) {
                        for (int src : graph.getIncoming(source_node)) {
                            if (kVerboseDbg) cout << "repetidos " << src << " -> " << source_node << endl;
                            graph.insertEdge(src, source_node + 1, 0);
                            graph.insertIncoming(src, source_node + 1);                  
                            graph.insertOutComing(src, source_node + 1);
                        }
                    }
                } else {
                    graph.insertEdge(target_node, source_node, 0);
                }
            }
        }
    }


    graph.markInitials(1);
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
    ifstream meuArquivo(nomeArquivo);
    if (!meuArquivo) {
        throw runtime_error("Arquivo " + nomeArquivo + " de kmers não encontrado");
    }

    if (detalhes) {
        cout << "Criando o grafo de De Bruijn" << endl;
    }

    string linha;
    size_t qtd = 0;

    while (getline(meuArquivo, linha)) {
        getline(meuArquivo, linha);
        if (detalhes) {
            cout << "Adicionando " << linha << " " << k << endl;
        }
        insertSequence(linha);
        qtd += linha.size() - k;
    }

    if (detalhes) {
        cout << "De Bruijn criado." << endl;
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
                dbg_gab.insertKmer(itr->first, 0);
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
    dbg_aux->insertKmer(kmer, 0);
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
                dbg_gab.insertKmer(itr->first, 0);
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
            if (kVerboseDbg) cout << "achei " << itr->second.node_in_sequence_graph  << " " << index << endl;
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
                        dbg.insertKmer(kmer, 0);
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