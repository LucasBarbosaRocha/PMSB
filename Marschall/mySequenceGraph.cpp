/*
    Classe mySequenceGraph
    Autor: Lucas B. Rocha
    Orientador(es): Said Sadique e Francisco Elói
    Ano: 2021
    FACOM: Doutorado em Ciência da Computação    
    Implementação: 
        construção de um grafo de sequências simples com lista de adjacências
*/
#include <bits/stdc++.h>
#define INF INT_MAX
using namespace std;

class SequenceGraph
{
private:
	int V; /* número de vértices */
	int k; /* k-mer do SequenceGraph De Bruijn */
	vector<pair<int, int> > *adj; /* Representacao do grafo */
	map<int, bool> visited;
	vector<int> iniciais;
	bool val;
    vector<string> *bases;

public:
	/* Construtor do grafo de sequências. Recebe como entrada a qtd. de vertices
	   e o comprimento do k-mer e devolve um grafo de sequências simples vazio.*/
	SequenceGraph(int V, int k); 

	/* funcao insere um vertice e um caractere no grafo de sequências simples */
	void insertNode(int v1, string base);

	/* funcao insere uma aresta com peso no grafo de sequências simples */
    void insertEdge(int u, int v, int wt);

	void alterarPesoAresta(int u, int v, int wt);

	/* funcao verifica se v2 é vizinho de v1  */
	bool isThereNeighbor(int v1, int v2);

	/* funcao imprime um grafo de sequências */
    void printGraph();

	/* funcao dfs recebe um vertice e um comprimento l e busca
	   um caminho de comprimento l no grafo  */
	void dfs(int v, int length);

	/* funcao recebe um vertice e um comprimento l e verifica
	   se existe um caminho de comprimento l no grafo  */
	bool isExistsPathLength(int v, int length);

	/* funcao verifica e marcas os vertices iniciais no grafo para execucao do algoritmo de Marschall  */
	void markInitials();

	/* funcao recebe um vertice v e altera se o vertice eh ou nao inicial */
	void alterarValorVerticeInicial(int v, int value);

	/* funcao recebe um vertice v e verifica se v eh inicial */
	int isInicial(int v);

	/* funcao recebe um vertice v e devolve o grau de saida */
	int getOutDegree(int v);

	/* funcao devolve a quantidade de vertice do grafo  */
	int getV();

	/* funcao devolve o k (comprimento do k-mer) */
	int getK();

	/* funcao recebe um vertice i e devolve o comeco da lista de adjacencia do vertice */
	vector<pair<int, int>>::iterator getAdjBegin(int i);

	/* funcao recebe um vertice i e devolve o final da lista de adjacencia do vertice */
	vector<pair<int, int>>::iterator getAdjEnd(int i);

	/* funcao recebe um vertice i e devolve o caractere associado ao vertice (o seu rótulo) */
	string getBase(int i);

	int verifyLabelExistisAndReturnVerticeIndex(string label);

	pair<vector<pair<int,string>>, int> dijkstra(int orig, int dest);


};

SequenceGraph::SequenceGraph(int V, int k)
{
	this->k = k;
	this->V = V; // atribui o número de vértices
	bases = new vector<string>[V];
    this->adj = new vector<pair<int,int>>[V];	
	for(int i = 0; i < V; i++)
		iniciais.push_back(0);
}

void SequenceGraph::insertNode(int v1, string base)
{
	// adiciona vértice v2 à lista de vértices adjacentes de v1
    bases[v1].push_back(base);
}

// To add an edge
void SequenceGraph::insertEdge(int u, int v, int wt)
{
	this->adj[u].push_back(make_pair(v, wt));
}

// To add an edge
void SequenceGraph::alterarPesoAresta(int u, int v, int wt)
{
	for (auto it = this->adj[u].begin(); it != this->adj[u].end(); it++)
	{
		if ((*it).first == u)
			(*it).second = wt;
	}
}

bool SequenceGraph::isThereNeighbor(int v1, int v2)
{
	for (auto it = adj[v1].begin(); it != adj[v1].end(); it++)
		if ((*it).first == v2)
			return true;
	return false;
}

// Print adjacency list representation ot graph
void SequenceGraph::printGraph()
{
	cout << "Sequence graph:" << endl;
	cout << "Qty. nodes: " << this->V << endl;
    for (int i = 0; i < this->V; i++)
    {
        cout << "(" <<  iniciais[i] << ")" << i << "-" << bases[i].front() << ": ";
		//cout << i << "-" << bases[i].front() << ": ";
        for (auto it = adj[i].begin(); it != adj[i].end(); it++)
            cout << (*it).first << "(" << (*it).second << ") ";
        cout << endl;
    }
}

void SequenceGraph::dfs(int v, int length)
{
    // Mark the current node as visited and
    // print it
    visited[v] = true;
    // cout << v << " ";

	if (this->k - 1 - length == 0)
		this->val = true;
	
    // Recur for all the vertices adjacent
    // to this vertex
    for (auto i = adj[v].begin(); i != adj[v].end(); ++i)
        if (!visited[(*i).first])
            dfs((*i).first, length - 1);
}

bool SequenceGraph::isExistsPathLength(int v, int length)
{
	this->val = false;
	dfs(v, length);
	this->visited.clear();
	if (this->val)
		return true;
	return false;
}
 
void SequenceGraph::markInitials()
{
    list<int>::iterator it;
    for (int i = 0; i < this->V; i++)
    {
        if(isExistsPathLength(i, this->k+1))
		{
			//ini[i] = true;
			this->alterarValorVerticeInicial(i, 1);           
		}
    }
}

int SequenceGraph::getOutDegree(int v)
{
	// basta retornar o tamanho da lista que é a quantidade de vizinhos
	return adj[v].size();
}

int SequenceGraph::isInicial(int v)
{
	if (this->iniciais[v] == 1)
		return 1;
	return 0;
}

void SequenceGraph::alterarValorVerticeInicial(int v, int value)
{
	this->iniciais[v] = value;
}

int SequenceGraph::getV()
{
	return this->V;
}

int SequenceGraph::getK()
{
	return this->k;
}

vector<pair<int, int>>::iterator SequenceGraph::getAdjBegin(int i)
{
	return this->adj[i].begin();
}

vector<pair<int, int>>::iterator SequenceGraph::getAdjEnd(int i)
{
	return this->adj[i].end();
}

string SequenceGraph::getBase(int i)
{
	return this->bases[i].front();
}

// Print adjacency list representation ot graph
int SequenceGraph::verifyLabelExistisAndReturnVerticeIndex(string label)
{
    for (int i = 0; i < this->V; i++)
    {
		if (bases[i].front().compare(label) == 0)
			return i;
    }
	return -1;
}

// Dijkstra
pair<vector<pair<int,string>>, int> SequenceGraph::dijkstra(int orig, int dest)
{

    // vetor de distâncias
    int V = this->V;

    int dist[V];
    int prev[V];
    /*
        vetor de visitados serve para caso o vértice já tenha sido
        expandido (visitado), não expandir mais
    */
    int visitados[V];

    // fila de prioridades de pair (distancia, vértice)
    priority_queue <pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;

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
			if (this->adj[u].size() > 0)
			{
				for(auto it = this->getAdjBegin(u); it !=  this->getAdjEnd(u); it++)
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

    list<string> induced_sequence;
    vector<pair<int,string>> saida;
    induced_sequence.push_back(this->getBase(dest));
    saida.push_back(make_pair(dest,this->getBase(dest)));

    for (int j = dest; j >= 0; j = prev[j])
    {
		if (prev[j] != -1)
		{
			induced_sequence.push_back(this->getBase(prev[j]));
			saida.push_back(make_pair(prev[j],this->getBase(prev[j])));
		}
    }

    // defvolve a distância mínima até o destino
    return make_pair(saida, dist[dest]);
}

// Driver code
/* int main()
{
	int V = 5;
    SequenceGraph SequenceGraph(V,3);
	SequenceGraph.insertEdge(0, 1, 10);
	SequenceGraph.insertEdge(0, 4, 20);
	SequenceGraph.insertEdge(1, 2, 30);
	SequenceGraph.insertEdge(1, 3, 40);
	SequenceGraph.insertEdge(2, 3, 40);
	SequenceGraph.insertEdge(3, 4, 40);
	//SequenceGraph.insertEdge(4, 0, 40);
	SequenceGraph.printGraph();

	int v1 = 2, v2 = 4, c = 1;
	cout << v1 << " -> " << v2 << "? ";
	if (SequenceGraph.isThereNeighbor(v1, v2))
		cout << "sim\n";
	else
		cout << "não\n";

	return 0;
} */
