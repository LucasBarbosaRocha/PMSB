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

private:
    int initialNode, endNode;

public:
    /* construtor da classe */
    Marschall(){};
    
    /* a função recebe dois caracteres A e B e os comparam, 
       devolve sub caso A != B, e 0 caso contrário */
    int w_sub(string caractere_grafo, string caractere_sequence);

    /* a função recebe um grafo de sequências simples G e uma sequência s,
       devolve um grafo de multicamadas com pesos nas arestas */
    SequenceGraph buildMultilayerGraph(SequenceGraph grafo, string sequence);

    /* a função de dijkstra recebe um grafo e dois vertices de origem e destino
       devolve a sequencia induzida pelo caminho mínimo e seu custo */
    pair<vector<pair<int,string>>, int>  dijkstra(SequenceGraph grafo, int orig, int dest);

    /* a função recebe dois vertices e salva o vertice inicial s 
       e vertice final t para a execução do dijkstra */
    void insertInitialAndEndNode(int v_initial, int v_end);

    /* devolve o vertice inicial para rodar no dijkstra */
    int getInitialNode();

    /* devolve o vertice final para rodar no dijkstra */
    int getEndNode();

    void shortestPath(SequenceGraph grafo, int src, int dest, int W);

    string verificaAresta(int u, int v, int tamGraph);

    pair<list<string>, string> showTraditionalMapping(vector<pair<int,string>> retorno, Hash deBruijnGraph, SequenceGraph traditionalGraph);
    pair<list<string>, string> showSimplifiedMapping(vector<pair<int,string>> retorno,  Hash deBruijnGraph, SequenceGraph simplifiedGraph);


};

int Marschall::w_sub(string caractere_grafo, string caractere_sequence)
{
    if (caractere_grafo.compare(caractere_sequence) == 0)
        return 0;
    return sub;
}

SequenceGraph Marschall::buildMultilayerGraph(SequenceGraph grafo, string sequence)
{
    int V = grafo.getV(), m = sequence.length(), vertice_atual = 0, vertice_inicial = 0, vertice_final, vertice_atual_aux, controle;
    int m_v = m * (V + 1) + 2; // quantidade de vertice do grafo multicamadas
    SequenceGraph m_grafo(m_v, grafo.getK());
    int *mapeamento;

    mapeamento = new (nothrow) int[V];
    this->sequenceGraphAndMulticamada = new (nothrow) vector<int>[m_v];  
    if (mapeamento == nullptr || this->sequenceGraphAndMulticamada == nullptr)
    {
        cerr << "error allocation multlayer graph" << endl;
        return m_grafo;
    }

    for (int i = 0; i <= m; i++)
    {      
        if (i == 0) // camada inicial
        {
            m_grafo.insertNode(vertice_atual, "s");
            vertice_atual++;
        } else {
            vertice_atual_aux = vertice_atual; 
            // dummy
            sequenceGraphAndMulticamada[vertice_atual_aux].push_back(-1);
            m_grafo.insertNode(vertice_atual_aux, "d");
            vertice_atual_aux++;
            // vertices
            for (int j = 0; j < V; j++)
            {
                string base = grafo.getBase(j);
                m_grafo.insertNode(vertice_atual_aux, base);
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
                    m_grafo.insertEdge(mapeamento[j], mapeamento[(*it).first], ins);
                }
            }

            if (i - 1 == 0)
            {
                // dummy
                m_grafo.insertEdge(vertice_inicial, vertice_atual, del);
                for (int j = 0; j < V; j++)
                {
                    // substituicao
                    if (grafo.isInicial(j))
                    {
                        m_grafo.insertEdge(vertice_inicial, mapeamento[j], w_sub(grafo.getBase(j), sequence.substr(i-1,1)));    
                    }              
                }   
            } else {
                int vertice_atual_camada_anterior = vertice_atual - (V + 1);
                int dif;

                // dummy esta em vertice_atual
                // delecao
                m_grafo.insertEdge(vertice_atual - (V + 1), vertice_atual, del);                
                for (int j = 0; j < V; j++)
                {
                    // substituicao
                    if (grafo.isInicial(j)) // j eh no grafo original
                        m_grafo.insertEdge(vertice_atual - (V + 1), mapeamento[j], w_sub(grafo.getBase(j), sequence.substr(i-1,1)));                  
                }                
           
                // outros
                for (int j = 0; j < V; j++)
                {
                    // delecao
                    m_grafo.insertEdge(mapeamento[j] - (V + 1), mapeamento[j], del);
                    // substituicao
                    for (auto it = grafo.getAdjBegin(j); it != grafo.getAdjEnd(j); it++)
                    {
                        m_grafo.insertEdge(mapeamento[j] - (V + 1), mapeamento[(*it).first], w_sub(grafo.getBase((*it).first), sequence.substr(i-1,1)));
                    }
                }                
            }
            vertice_atual = vertice_atual_aux; // atualizando o vertice atual
        }
    }

    // criar ultimo vertice
    vertice_final = vertice_atual;
    m_grafo.insertNode(vertice_final, "t");
    // dummy
    m_grafo.insertEdge(vertice_final - (V + 1), vertice_final, 0);
    // outros

    for (int j = 0; j < V; j++)
    {
        m_grafo.insertEdge(vertice_final - (V + 1) + (j+1), vertice_final, 0);
    }
    // criar grafo multicamadas

    this->insertInitialAndEndNode(vertice_inicial, vertice_final);

    // deletando vetores sem utilizacao
    delete [] mapeamento;

    return m_grafo;
}

// Dijkstra
pair<vector<pair<int,string>>, int> Marschall::dijkstra(SequenceGraph grafo, int orig, int dest)
{
    // vetor de distâncias
    int V = grafo.getV();
    list<string> induced_sequence;
    vector<pair<int,string>> saida;
    int *dist;
    int *prev;
    /*
        vetor de visitados serve para caso o vértice já tenha sido
        expandido (visitado), não expandir mais
    */
    int *visitados;
    // fila de prioridades de pair (distancia, vértice)
    priority_queue <pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;

    dist = new (nothrow) int[V];
    prev = new (nothrow) int[V];
    visitados = new (nothrow) int[V];
    //new (nothrow) list<int>::iterator[V];
  
    if (prev == nullptr || dist == nullptr || visitados == nullptr)
    {
        cerr << "error allocation" << endl;
        return make_pair(saida, -1);
    }

    // inicia o vetor de distâncias e visitados
    for(int i = 0; i < V; i++)
    {
        dist[i] = INF;
        prev[i] = -1;
        visitados[i] = false;

    }

    // a distância de orig para orig é 0
    dist[orig] = 0;

    // insere na fila
    pq.push(make_pair(dist[orig], orig));

    // loop do algoritmo
    while(!pq.empty())
    {
        pair<int, int> p = pq.top(); // extrai o pair do topo
        int u = p.second; // obtém o vértice do pair
        pq.pop(); // remove da fila

        // verifica se o vértice não foi expandido
        if(visitados[u] == false)
        {
            // marca como visitado
            visitados[u] = true;
            // percorre os vértices "v" adjacentes de "u" se existire size() > 0
			if (grafo.getOutDegree(u) > 0)
			{
                // percorre os vértices "v" adjacentes de "u"
                for(auto it = grafo.getAdjBegin(u); it !=  grafo.getAdjEnd(u); it++)
                {
                    // obtém o vértice adjacente e o custo da aresta
                    int v = it->first;
                    int custo_aresta = it->second;

                    // relaxamento (u, v)
                    if(dist[v] > (dist[u] + custo_aresta))
                    {
                        // atualiza a distância de "v" e insere na fila
                        dist[v] = dist[u] + custo_aresta;
                        prev[v] = u;
                        pq.push(make_pair(dist[v], v));
                    }
                }
            }
        }
    }

    induced_sequence.push_back(grafo.getBase(dest));
    saida.push_back(make_pair(dest,grafo.getBase(dest)));
    for (int j = dest; j > 0; j = prev[j])
    {
        if (prev[j] != -1)
        {
            induced_sequence.push_back(grafo.getBase(prev[j]));
            saida.push_back(make_pair(prev[j],grafo.getBase(prev[j])));
        }
    }

    /*for (auto it = saida.begin(); it != saida.end(); it++)
    {
        cout << (*it).first << ":" << (*it).second << " ";
    }
    cout << endl;*/
    // retorna a distância mínima até o destino

    int custo = dist[dest];

        delete [] visitados;
        delete [] prev;
        delete [] dist;
    

    return make_pair(saida, custo);
}

void Marschall::insertInitialAndEndNode(int v_initial, int v_end)
{
    this->initialNode = v_initial;
    this->endNode = v_end;
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
    cout << dist2[dest] << endl;

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

pair<list<string>, string> Marschall::showTraditionalMapping(vector<pair<int,string>> retorno, Hash deBruijnGraph, SequenceGraph traditionalGraph)
{
    int primeiro = 0, indice, anterior = 0, details = 0, k = traditionalGraph.getK(), kmer_count = 0;
    string aux, tmp, baseAnterior, kmer_aux = "";
    list<string> kmers;

    //cout << "mapeando " << endl;

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

            tmp = this->verificaAresta((*it).first, anterior, traditionalGraph.getV());
            if (details == 1)
                cout << "(" << tmp << ") ";
            anterior = (*it).first;    
            indice = this->sequenceGraphAndMulticamada[(*it).first].front(); 

            if (indice != -1)
            {   
                // auto kmer = kmerAndNode.at(indice);      
                auto kmer = deBruijnGraph.findKmerByIndex(indice);
                if (details == 1)
                    cout << (*it).second << "(" << kmer << ") <-";

                if (kmer.compare(kmer_aux) != 0 || kmer_count == 0)
                {
                    kmers.push_front(kmer);
                    kmer_count = 0;
                    kmer_aux = kmer;
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

    /*for (auto a = kmers.begin(); a != kmers.end(); a++)
    {
        cout << (*a) << " ";
    }
    cout << endl;*/

    return  make_pair(kmers,aux.substr(0, aux.length() - 1));
}

pair<list<string>, string> Marschall::showSimplifiedMapping(vector<pair<int,string>> retorno, Hash deBruijnGraph, SequenceGraph simplifiedGraph)
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

