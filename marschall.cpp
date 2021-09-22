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

using namespace std;

class Marschall
{
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
    pair<list<string>, int>  dijkstra(SequenceGraph grafo, int orig, int dest);

    /* a função recebe dois vertices e salva o vertice inicial s 
       e vertice final t para a execução do dijkstra */
    void insertInitialAndEndNode(int v_initial, int v_end);

    /* devolve o vertice inicial para rodar no dijkstra */
    int getInitialNode();

    /* devolve o vertice final para rodar no dijkstra */
    int getEndNode();
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
    int mapeamento[V];

    for (int i = 0; i <= m; i++)
    {      
        if (i == 0) // camada inicial
        {
            m_grafo.insertNode(vertice_atual, "s");
            vertice_atual++;
        } else {
            vertice_atual_aux = vertice_atual; 
            // dummy
            m_grafo.insertNode(vertice_atual_aux, "d");
            vertice_atual_aux++;
            // vertices
            for (int j = 0; j < V; j++)
            {
                string base = grafo.getBase(j);
                m_grafo.insertNode(vertice_atual_aux, base);
                mapeamento[j] = vertice_atual_aux;
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
    return m_grafo;
}

// Dijkstra
pair<list<string>, int> Marschall::dijkstra(SequenceGraph grafo, int orig, int dest)
{
    // vetor de distâncias
    int V = grafo.getV();
    int dist[V];
    int prev[V];

    /*
        vetor de visitados serve para caso o vértice já tenha sido
        expandido (visitado), não expandir mais
    */
    int visitados[V];

    // fila de prioridades de pair (distancia, vértice)
    priority_queue < pair<int, int>,
                    vector<pair<int, int> >, greater<pair<int, int> > > pq;

    // inicia o vetor de distâncias e visitados
    for(int i = 0; i < V; i++)
    {
        dist[i] = INF;
        prev[-1];
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

    list<string> induced_sequence;
    induced_sequence.push_back(grafo.getBase(dest));
    for (int j = dest; j > 0; j = prev[j])
    {
        induced_sequence.push_back(grafo.getBase(prev[j]));
    }

    // retorna a distância mínima até o destino
    return make_pair(induced_sequence, dist[dest]);
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

int main()
{
    string sequence = "CGA";
    int k = 3;

    // insert the kmers into the hash table
    Hash h(k); // 7 is count of buckets in
                // hash table
    
    Marschall m;
    
    string nomeArquivo = "kmers4.txt";
    
    h.populateGraph(nomeArquivo, false); 
    auto grafo = h.dbgToSequenceGraph_1();
    auto m_grafo = m.buildMultilayerGraph(grafo, sequence);
    auto retorno = m.dijkstra(m_grafo, m.getInitialNode(), m.getEndNode());

    // m_grafo.printGraph();

    for (auto it = retorno.first.begin(); it != retorno.first.end(); it++)
    {
        cout << (*it) << " <- ";
    }
    cout << endl;
    cout << "Cost " << retorno.second << endl; 

    auto grafo2 = h.dbgToSequenceGraph_2();
    auto m_grafo2 = m.buildMultilayerGraph(grafo2, sequence);
    retorno = m.dijkstra(m_grafo2,  m.getInitialNode(), m.getEndNode());
    for (auto it = retorno.first.begin(); it != retorno.first.end(); it++)
    {
        cout << (*it) << " <- ";
    }
    cout << endl;
    cout << "Cost " << retorno.second << endl; 

return 0;
}