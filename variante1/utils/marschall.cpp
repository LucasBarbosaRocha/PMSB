/*
    Classe Marschall
    Autor: Lucas B. Rocha
    Orientador(es): Said Sadique e Francisco Elói
    Ano: 2021
    FACOM: Doutorado em Ciência da Computação    
    Implementação: 
        construção do grafo de multicamadas e implementação do dijkstra. 
        Podemos executar o dijkstra e determinar o menor caminho no grafo de multicamadas
        conforma a ideia de Rautiainen and Marschall
    Artigo: Aligning sequences to general graphs in O(V + mE) time
    Autores: Rautiainen and Marschall
*/
#include <iostream>
#include "myDBGgraph.cpp"
#include <queue>
#include <bits/stdc++.h>

#define sub 1
#define ins 1
#define del 1
#define INF INT_MAX
string sequence;
string nameArchive;
int k;

using namespace std;

class Marschall
{
public:
    vector<int> *sequenceGraphAndMulticamada;
    SequenceGraph m_sequenceGraph;

private:
    int initialNode, endNode;

public:
    /* construtor da classe */
    Marschall() : sequenceGraphAndMulticamada(nullptr) {};
    ~Marschall() { delete[] sequenceGraphAndMulticamada; };
    
    /* a função recebe dois caracteres A e B e os comparam, 
       devolve sub caso A != B, e 0 caso contrário */
    int w_sub(string caractere_grafo, string caractere_sequence);

    /* a função recebe um grafo de sequências simples G e uma sequência s,
       devolve um grafo de multicamadas com pesos nas arestas */
    void buildMultilayerGraph(SequenceGraph &grafo, const string &sequence);

    /* a função de dijkstra recebe um grafo e dois vertices de origem e destino
       devolve a sequencia induzida pelo caminho mínimo e seu custo */
    pair<vector<pair<int,string>>, int> dijkstra(SequenceGraph &grafo, int orig, int dest, int limite);

    /* a função recebe dois vertices e salva o vertice inicial s
       e vertice final t para a execução do dijkstra */
    void insertInitialAndEndNode(int v_initial, int v_end);

    void invertInitialAndEndNode();

    /* devolve o vertice inicial para rodar no dijkstra */
    int getInitialNode();

    /* devolve o vertice final para rodar no dijkstra */
    int getEndNode();

    void shortestPath(SequenceGraph grafo, int src, int dest, int W);

    string verificaAresta(int u, int v, int tamGraph);

    pair<list<string>, string> showTraditionalMapping(const vector<pair<int,string>> &retorno, Hash &deBruijnGraph, SequenceGraph &traditionalGraph);
    pair<list<string>, string> showSimplifiedMapping(const vector<pair<int,string>> &retorno, Hash &deBruijnGraph, SequenceGraph &simplifiedGraph);


};

int Marschall::w_sub(string caractere_grafo, string caractere_sequence)
{
    if (caractere_grafo.compare(caractere_sequence) == 0)
        return 0;
    return sub;
}

void Marschall::buildMultilayerGraph(SequenceGraph &grafo, const string &sequence)
{
    int V = grafo.getV(), m = sequence.length(), vertice_atual = 0, vertice_inicial = 0, vertice_final, vertice_atual_aux, controle;
    int m_v = m * (V + 1) + 2; // quantidade de vertice do grafo multicamadas
    //SequenceGraph m_grafo(m_v, grafo.getK());
    m_sequenceGraph.initilizeSequenceGraph(m_v, grafo.getK());
    int *mapeamento;

    mapeamento = new (nothrow) int[V];
    delete[] this->sequenceGraphAndMulticamada;
    this->sequenceGraphAndMulticamada = new (nothrow) vector<int>[m_v];
    if (mapeamento == nullptr || this->sequenceGraphAndMulticamada == nullptr)
    {
        cerr << "error allocation multlayer graph" << endl;
    }

    for (int i = 0; i <= m; i++)
    {      
        if (i == 0) // camada inicial
        {
            m_sequenceGraph.insertNode(vertice_atual, "s");
            vertice_atual++;
        } else {
            vertice_atual_aux = vertice_atual; 
            // dummy
            sequenceGraphAndMulticamada[vertice_atual_aux].push_back(-1);
            m_sequenceGraph.insertNode(vertice_atual_aux, "d");
            vertice_atual_aux++;
            // vertices
            for (int j = 0; j < V; j++)
            {
                string base = grafo.getBase(j);
                m_sequenceGraph.insertNode(vertice_atual_aux, base);
                mapeamento[j] = vertice_atual_aux;
                sequenceGraphAndMulticamada[vertice_atual_aux].push_back(j);
                vertice_atual_aux++;               
            }
            
            // arestas adjacentes
            for (int j = 0; j < V; j++)
            {
                for (auto it = grafo.getAdjBegin(j); it != grafo.getAdjEnd(j); it++)
                {
                    // insercao
                    m_sequenceGraph.insertEdge(mapeamento[j], mapeamento[(*it).first], ins);
                }
            }

            if (i - 1 == 0)
            {
                // dummy
                m_sequenceGraph.insertEdge(vertice_inicial, vertice_atual, del);
                for (int j = 0; j < V; j++)
                {
                    // substituicao
                    if (grafo.isInicial(j))
                    {
                        //cout << "COISA compara " << vertice_inicial << " -> " << mapeamento[j] << " " << j << " " << grafo.getBase(j) << " =? " << sequence.substr(i-1,1) << " custo " << w_sub(grafo.getBase(j), sequence.substr(i-1,1)) << endl;
                        m_sequenceGraph.insertEdge(vertice_inicial, mapeamento[j], w_sub(grafo.getBase(j), sequence.substr(i-1,1)));    
                    }              
                }   
            } else {
                int vertice_atual_camada_anterior = vertice_atual - (V + 1);
                int dif;

                // dummy esta em vertice_atual
                // delecao
                m_sequenceGraph.insertEdge(vertice_atual - (V + 1), vertice_atual, del);                
                for (int j = 0; j < V; j++)
                {
                    // substituicao
                    if (grafo.isInicial(j)) // j eh no grafo original 
                    {
                        //cout << "COISA compara " << vertice_atual - (V + 1) << " -> " << mapeamento[j] << " " << j << " " << grafo.getBase(j) << " =? " << sequence.substr(i-1,1) << " custo " << w_sub(grafo.getBase(j), sequence.substr(i-1,1)) << endl;
                        m_sequenceGraph.insertEdge(vertice_atual - (V + 1), mapeamento[j], w_sub(grafo.getBase(j), sequence.substr(i-1,1)));                  
                    }
                }                
           
                // outros
                for (int j = 0; j < V; j++)
                {
                    // delecao
                    m_sequenceGraph.insertEdge(mapeamento[j] - (V + 1), mapeamento[j], del);
                    // substituicao
                    for (auto it = grafo.getAdjBegin(j); it != grafo.getAdjEnd(j); it++)
                    {
                        m_sequenceGraph.insertEdge(mapeamento[j] - (V + 1), mapeamento[(*it).first], w_sub(grafo.getBase((*it).first), sequence.substr(i-1,1)));
                    }
                }                
            }
            vertice_atual = vertice_atual_aux; // atualizando o vertice atual
        }
    }

    // criar ultimo vertice
    vertice_final = vertice_atual;
    m_sequenceGraph.insertNode(vertice_final, "t");
    // dummy
    m_sequenceGraph.insertEdge(vertice_final - (V + 1), vertice_final, 0);
    // outros

    for (int j = 0; j < V; j++)
    {
        m_sequenceGraph.insertEdge(vertice_final - (V + 1) + (j+1), vertice_final, 0);
    }
    // criar grafo multicamadas
    this->insertInitialAndEndNode(vertice_inicial, vertice_final);

    // deletando vetores sem utilizacao
    delete [] mapeamento;
}

// Dijkstra
pair<vector<pair<int, string>>, int> Marschall::dijkstra(SequenceGraph &grafo, int orig, int dest, int limite = -1) {
    int V = grafo.getV();
    vector<pair<int, string>> saida;
    vector<int> dist(V, INF), prev(V, -1), visitados(V, false);

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    dist[orig] = 0;
    pq.emplace(dist[orig], orig);

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        if (!visitados[u]) {
            visitados[u] = true;
            if (grafo.getOutDegree(u) > 0) {
                for (auto it = grafo.getAdjBegin(u); it != grafo.getAdjEnd(u); ++it) {
                    int v = it->first;
                    int custo_aresta = it->second;
                    if (dist[v] > dist[u] + custo_aresta) {
                        dist[v] = dist[u] + custo_aresta;
                        prev[v] = u;
                        pq.emplace(dist[v], v);
                    }
                }
            }
        }
    }

    if (limite > -1) {
        cout << "Custo\n" << dist[orig] << "\n";
        for (int i = 1; i < V - limite; i += limite) {
            for (int j = 0; j < limite; ++j) {
                cout << dist[i + j] << " ";
            }
            cout << "\n";
        }
        cout << dist[dest] << "\n\n";
    }

    list<string> induced_sequence;
    induced_sequence.push_back(grafo.getBase(dest));
    saida.push_back(make_pair(dest,grafo.getBase(dest)));
    for (int j = dest; j > 0; j = prev[j]) {
        int prevIndex = prev[j];
        if (prevIndex != -1) {
            auto base = grafo.getBase(prevIndex);
            induced_sequence.push_back(base);
            saida.emplace_back(prevIndex, base);
        }
    }
    
    return make_pair(saida, dist[dest]);
}

void Marschall::insertInitialAndEndNode(int v_initial, int v_end)
{
    this->initialNode = v_initial;
    this->endNode = v_end;
}

void Marschall::invertInitialAndEndNode()
{
    auto aux = this->initialNode;
    this->initialNode = this->endNode;
    this->endNode = aux;
}

int Marschall::getInitialNode()
{
    return this->initialNode;
}

int Marschall::getEndNode()
{
    return this->endNode;
}

// Prints shortest paths from src to all other vertices.
// W is the maximum weight of an edge
void Marschall::shortestPath(SequenceGraph grafo, int src, int dest, int W)
{
    /* With each distance, iterator to that vertex in
       its bucket is stored so that vertex can be deleted
       in O(1) at time of updation. So
    dist[i].first = distance of ith vertex from src vertex
    dits[i].second = iterator to vertex i in bucket number */
    int V = grafo.getV();

    //vector<pair<int, list<int>::iterator>> dist(V);
    //int prev[V];

    int *prev;
    //vector<pair<int, list<int>::iterator>> *dist;
    int *dist2;
    list<int>::iterator *buckets;


    prev = new (nothrow) int[V];
    dist2 = new (nothrow) int[V];
    buckets = new (nothrow) list<int>::iterator[V];
  
    if (prev == nullptr || dist2 == nullptr || buckets == nullptr)
    {
        cerr << "error allocation prev vector" << endl;
        return;
    }

    //vector<pair<int, list<int>::iterator>> dist(V);
    // Initialize all distances as infinite (INF)
    for (int i = 0; i < V; i++)
    {
        dist2[i] = INF;
        prev[i] = -1;
    }

    // Create buckets B[].
    // B[i] keep vertex of distance label i
    //list<int> B[W * V + 1];
    unordered_map<int, list<int>> B;  
    B[0].push_back(src);
    dist2[src] = 0;
    
    int idx = 0;
    while (1)
    {
        // Go sequentially through buckets till one non-empty
        // bucket is found
        while (B[idx].size() == 0 && idx < W*V)
            idx++;
  
        // If all buckets are empty, we are done.
        if (idx == W * V)
            break;
  
        // Take top vertex from bucket and pop it
        int u = B[idx].front();
        B[idx].pop_front();
        // Process all adjacents of extracted vertex 'u' and
        // update their distanced if required.
        for (auto i = grafo.getAdjBegin(u); i != grafo.getAdjEnd(u); ++i)
        {
            int v = (*i).first;
            int weight = (*i).second;
  
            int du = dist2[u];
            int dv = dist2[v];
  
            // If there is shorted path to v through u.
            if (dv > du + weight)
            {
                // If dv is not INF then it must be in B[dv]
                // bucket, so erase its entry using iterator
                // in O(1)
                if (dv != INF)
                    B[dv].erase(buckets[v]);
  
                //  updating the distance
                dist2[v] = du + weight;
                dv = dist2[v];
                prev[v] = u;
  
                // pushing vertex v into updated distance's bucket
                B[dv].push_front(v);
  
                // storing updated iterator in dist[v].second
                buckets[v] = B[dv].begin();
            }
        }
    }  
    // Print shortest distances stored in dist[]
    //printf("Vertex   Distance from Source\n");
    //for (int i = 0; i < V; ++i)
        //printf("%d     %d\n", i, dist[i].first);

    list<string> induced_sequence;
    vector<pair<int,string>> saida;
    induced_sequence.push_back(grafo.getBase(dest));
    saida.push_back(make_pair(dest,grafo.getBase(dest)));
    for (int j = dest; j > 0; j = prev[j])
    {
        induced_sequence.push_back(grafo.getBase(prev[j]));
        saida.push_back(make_pair(prev[j],grafo.getBase(prev[j])));
    }

    /*for (auto it = saida.begin(); it != saida.end(); it++)
    {
        cout << (*it).first << ":" << (*it).second << " ";
    }
    cout << endl; */

    /* Liberando memória */
    B.clear();
    delete [] dist2;
    delete [] buckets;
    delete [] prev;
}

int verificaEntrada(int argc, char *argv[])
{
    string aux; 
    if (argc == 1)
    {
        cout << "Error: digite -help" << endl;
        return 1;
    }

    if (argc == 2)
    {
        aux = argv[1];
        if (aux.compare("-help") == 0)
            cout << "-s sequence -g graph -k int" << endl;
        return 1;
    }

    if (argc == 7)
    {
        aux = argv[1];
        if (aux.compare("-s") == 0)
            sequence = argv[2];
        else
        {
            cout << "Error: digite -help" << endl;
            return 1;
        }
        aux = argv[3];
        if (aux.compare("-g") == 0)
            nameArchive = argv[4];
        else
        {
            cout << "Error: digite -help" << endl;
            return 1;
        }
        aux = argv[5];
        if (aux.compare("-k") == 0)
            k = atoi(argv[6]);
        else
        {
            cout << "Error: digite -help" << endl;
            return 1;
        }
        //cout << sequence << " " << nameArchive << endl;
        return 0;
    }

    cout << "Error: digite -help" << endl;
    return 1;
}

string Marschall::verificaAresta(int u, int v, int tamGraph)
{
    int lim = u + (tamGraph + 1);
    //cout << u << " -> " << v << ":" << tamGraph << "lim: " << lim << " "; 

    if (v == lim)
        return "del";
    if (v >= lim - ((tamGraph+1)/2) && v < lim)
        return "sub";
    if (v > lim)
        return "sub";
    return "ins";
}

pair<list<string>, string> Marschall::showTraditionalMapping(const vector<pair<int, string>> &retorno, Hash &deBruijnGraph, SequenceGraph &traditionalGraph) {
    int anterior = 0, k = traditionalGraph.getK(), kmer_count = 0;
    string aux, tmp, baseAnterior, kmer_aux = "";
    list<string> kmers;

    for (auto it = retorno.begin(); it != retorno.end(); ++it) {
        tmp = "";
        
        if (it == retorno.begin()) {
            anterior = it->first;
            baseAnterior = it->second;
            continue;
        } 

        if (it == retorno.end() - 1) {
            tmp = anterior == 1 ? "del" : "sub";
            aux = baseAnterior + aux;
        } else {
            tmp = this->verificaAresta(it->first, anterior, traditionalGraph.getV());
            anterior = it->first;    
            int indice = this->sequenceGraphAndMulticamada[it->first].front(); 

            if (indice != -1) {   
                auto kmer = deBruijnGraph.findKmerByIndex(indice);
                
                if (kmer.compare(kmer_aux) != 0 || kmer_count == 0) {
                    kmers.push_front(kmer);
                    kmer_count = 0;
                    kmer_aux = kmer;
                }
                kmer_count++;
            }

            aux = (tmp == "del") ? "-" + aux : baseAnterior + aux;
            baseAnterior = (indice == -1) ? "-" : it->second;
        } 
    }

    return make_pair(kmers, aux.substr(0, aux.length() - 1));
}

pair<list<string>, string> Marschall::showSimplifiedMapping(const vector<pair<int,string>> &retorno, Hash &deBruijnGraph, SequenceGraph &simplifiedGraph)
{
    int primeiro = 0, indice, anterior = 0, details = 0, k = simplifiedGraph.getK(), kmer_count = 0;
    string aux, tmp, baseAnterior, kmer_aux = "", kmerMapeado;
    list<string> kmers;

    if (details == 1)
    {
        for (auto it = retorno.begin(); it != retorno.end(); it++)
        {
            cout << (*it).second << " <- ";
            aux = (*it).second + aux;
        }
        cout << endl;
        cout << aux << endl; 
    }

    for (auto it = retorno.begin(); it != retorno.end(); it++)
    {
        tmp = "";
        if (it == retorno.end() - 1)
        {
            if (anterior == 1)
                tmp = "del";
            else
                tmp = "sub";
            aux = baseAnterior + aux;
            if (details == 1)
                cout << "(" << tmp << ") ";
        }else if (it == retorno.begin())
        {
            anterior = (*it).first;
            baseAnterior = (*it).second;
        }
        else
        {

            tmp = this->verificaAresta((*it).first, anterior, simplifiedGraph.getV());
            if (details == 1)
                cout << "(" << tmp << ") ";
            anterior = (*it).first;    
            indice = this->sequenceGraphAndMulticamada[(*it).first].front(); 

            if (indice != -1)
            {   
                auto kmer = deBruijnGraph.findKmerByIndex(indice);       
                if (details == 1)
                    cout << (*it).second << "(" << kmer << ") <-";

                if (kmer.find("$") == 0)
                    kmerMapeado = deBruijnGraph.findKmerBySpecialKmer(kmer); // mapeando um kmer do G'_k no G_K
                else
                    kmerMapeado = kmer;

                if (kmerMapeado.compare(kmer_aux) != 0 || kmer_count == 0)
                {
                    kmers.push_front(kmerMapeado);
                    kmer_count = 0;
                    kmer_aux = kmerMapeado;
                }
                kmer_count++;
            }

            if (tmp == "del")
            {
                aux = "-" + aux;
            }else
            {
                aux = baseAnterior + aux;
                baseAnterior = (*it).second; 
            } 

            if (indice == -1)
                baseAnterior = "-";    
        } 
    }
    return  make_pair(kmers,aux.substr(0, aux.length() - 1));
}

