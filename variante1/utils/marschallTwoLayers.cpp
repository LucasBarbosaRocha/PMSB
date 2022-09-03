#include <iostream>
//#include "myDBGgraph.cpp"
#include "myBifrost/my_bifrost.cpp"
#include <queue>
#include <bits/stdc++.h>

#define sub 1
#define ins 1
#define del 1

using namespace std;

class MarschallTwoLayers
{
private:
    int w(SequenceGraph grafo, int layer_u, int node_u, string base_string_v, int layer_v, int node_v, string base_node_v);
    vector<pair<int,int>> algoritmo_3(SequenceGraph sequenceGraph, vector<pair<int,int>> currentLayer, unordered_map<int, int> currentLayerBkp);
    vector<pair<int,int>> algoritmo_2(SequenceGraph sequenceGraph, string sequence, int layer, vector<pair<int,int>> PreviousLayer, vector<pair<int,int>> CurrentLayer, unordered_map<int, int> CurrentLayerBkp);

public:
    /* construtor da classe */
    MarschallTwoLayers(){};
    int toCalculateTheCost(SequenceGraph sequenceGraph, string sequence);   
};

bool sortbysec(const pair<int,int> &a, const pair<int,int> &b)
{
    return (a.second < b.second);
}

int MarschallTwoLayers::w(SequenceGraph grafo, int layer_u, int node_u, string base_string_v, int layer_v, int node_v, string base_node_v)
{
    // camada 0
    if (layer_u == 0) 
    {
        if (base_string_v.compare(base_node_v))
            return sub;
        else
            return 0;        
    }

    // camada > 0
    // sub
    if (layer_u == layer_v - 1 && grafo.isThereNeighbor(node_u,node_v))
    {
        if(base_string_v.compare(base_node_v))
            return sub;
        else
            return 0;
    }

    // ins
    if (layer_u == layer_v && grafo.isThereNeighbor(node_u,node_v))
    {
            return ins;
    }

    // del
    if (layer_u == layer_v - 1 && !grafo.isThereNeighbor(node_u,node_v))
    {
            return ins;
    }

    return 1;
}

vector<pair<int,int>> MarschallTwoLayers::algoritmo_3(SequenceGraph sequenceGraph, vector<pair<int,int>> currentLayer, unordered_map<int, int> currentLayerBkp)
{
    int V = sequenceGraph.getV(), resolved[V], x;
    vector<pair<int, int>>::iterator it;
    queue<int> q1, q2;

    for (int i = 0; i < V; i++)
        resolved[i] = 0;

    for (int i = 0; i < V; i++)
    {
        q1.push(currentLayer[i].first);
    }
    
    while (!q1.empty() || !q2.empty())
    {
        if (q1.size() < 1)
        {
            x = q2.front();
            q2.pop();
        }else if (q2.size() < 1)
        {
            x = q1.front();
            q1.pop();
        }else if (q1.front() < q2.front())
        {
            x = q1.front();
            q1.pop();
        }
        else
        {
            x = q2.front();
            q2.pop();
        }
        int pos_x = currentLayerBkp[x];
        if(resolved[pos_x] == 0)
        {
            resolved[pos_x] = 1;
            for (it = sequenceGraph.getAdjBegin(x); it != sequenceGraph.getAdjEnd(x); it++)
            {
                int pos_it = currentLayerBkp[(*it).first];
                if(currentLayer[pos_it].second > currentLayer[pos_x].second + ins)
                {
                    currentLayer[pos_it].second = currentLayer[pos_x].second + ins;
                    q2.push((*it).first);
                }
            }
        }
    }
    
    return currentLayer;
}

vector<pair<int,int>> MarschallTwoLayers::algoritmo_2(SequenceGraph sequenceGraph, string sequence, int layer, vector<pair<int,int>> PreviousLayer, vector<pair<int,int>> CurrentLayer, unordered_map<int, int> CurrentLayerBkp)
{
    int V = sequenceGraph.getV(), cost;
    vector<pair<int, int>>::iterator it;

    for (int p = 0; p < V; p++)
    {
        int j = PreviousLayer[p].first; // j is first node ordenado e p is the position of j in the vector

        // vizinhos do grafo
        cost = 0;
        for (it = sequenceGraph.getAdjBegin(j); it != sequenceGraph.getAdjEnd(j); it++)
        {
            // (pos e *it) e (p e j)  
            int pos = CurrentLayerBkp[(*it).first]; // *it eh node adj e vamos pegar a posicao desse noh na camada atual
            if(sequenceGraph.isThereNeighbor(j, (*it).first))
            {
                cost = w(sequenceGraph, layer, j, sequenceGraph.getBase((*it).first), layer+1, (*it).first, sequence.substr(layer, 1));
                if(CurrentLayer[pos].second >= PreviousLayer[p].second + cost)
                {
                    CurrentLayer[pos].second = PreviousLayer[p].second + cost;
                }
            }
        }
        // vizinho imediato na camada debaixo
        cost = w(sequenceGraph, layer, j, sequenceGraph.getBase(j), layer+1, j, sequence.substr(layer, 1));
        int pos = CurrentLayerBkp[j]; // *it eh node adj e vamos pegar a posicao desse noh na camada atual
        if(CurrentLayer[pos].second > PreviousLayer[p].second + cost)
        {
            CurrentLayer[pos].second = PreviousLayer[p].second + cost;
        }
    }
    sort(CurrentLayer.begin(), CurrentLayer.end(), sortbysec); 

    return CurrentLayer;
}

int MarschallTwoLayers::toCalculateTheCost(SequenceGraph sequenceGraph, string sequence)
{
    int V = sequenceGraph.getV(), m = sequence.length();
    vector<pair<int,int>> PreviousLayer, CurrentLayer; // node and cost
    unordered_map<int, int> CurrentLayerBkp; // node and position


    for (int i = 0; i < V; i++)
    {
        PreviousLayer.push_back(make_pair(i,0));
        CurrentLayerBkp[i] = i;
        CurrentLayer.push_back(make_pair(i,INT_MAX));
    }
    
    /* Começo Camada 0 */
    for (int i = 0; i < V; i++)
    {
        if(sequenceGraph.isInicial(i))
        {
            if(CurrentLayer[i].second > PreviousLayer[i].second + w(sequenceGraph, 0, i, sequenceGraph.getBase(i),1, 0, sequence.substr(0, 1)))
                CurrentLayer[i].second = PreviousLayer[i].second + w(sequenceGraph, 0, i, sequenceGraph.getBase(i),1, 0, sequence.substr(0, 1));
        }
    }
    sort(CurrentLayer.begin(), CurrentLayer.end(), sortbysec);

    for (int i = 0; i < V; i++)
    {
        PreviousLayer[i] = CurrentLayer[i]; 
        int pos = PreviousLayer[i].first;
        CurrentLayerBkp[pos] = i;
    }

    PreviousLayer = algoritmo_3(sequenceGraph, CurrentLayer, CurrentLayerBkp);
    /* Fim Camada 0 */

    for (int i = 1; i < m; i++)
    {
        for (int j = 0; j < V; j++)
            CurrentLayer[j].second = INT_MAX; 

        //algoritmo 2
        CurrentLayer = algoritmo_2(sequenceGraph, sequence, i, PreviousLayer, CurrentLayer, CurrentLayerBkp);
        // copiando as distancias das camadas e os indices
        for (int i = 0; i < V; i++)
        {
            PreviousLayer[i] = CurrentLayer[i]; 
            int pos = PreviousLayer[i].first;
            CurrentLayerBkp[pos] = i;
        }
        PreviousLayer = algoritmo_3(sequenceGraph, CurrentLayer, CurrentLayerBkp);
    }

    int menor = INT_MAX;
    for (int i = 0; i < V; i++)
    {
        // cout << PreviousLayer[i].second << " < " << menor << endl;
        if (PreviousLayer[i].second >= 0 && PreviousLayer[i].second < menor)
            menor = PreviousLayer[i].second;
    }
    return menor;
}

/*
int main(int argc, char *argv[])
{
    string line;
    if(utils.verifyData(argc, argv) == 1)
        exit (0); 
        
    Hash h(utils.k);   
    ifstream file(utils.nameSequenceArchive);
    h.populateGraph(utils.nameArchive, false);       
    h.dbgToSimplifiedSequenceGraph(0);
    while(getline(file, line))
    {
        //utils.readSequence(utils.nameSequenceArchive);  
        getline(file, line);
        cout << "Size L.Read " << line.size() << endl;
        utils.sequence = line;
        // mapeamento
        int retorno = toCalculateTheCost(h.sequenceGraph, utils.sequence);
        cout << retorno << endl << endl;
    }     
return 0;
}
*/