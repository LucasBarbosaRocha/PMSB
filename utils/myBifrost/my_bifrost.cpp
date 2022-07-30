//#include "bifrost/src/CompactedDBG.hpp"
#include "myData.cpp"
#include "../mySequenceGraph.cpp"

#define sub 1
#define ins 1
#define del 1

class my_Bifrost
{
private:
    int k;
    vector<int> *sequenceGraphAndMulticamada;
    CompactedDBG<myData> cdbg;
    int hammingDistance(string kmer1, string kmer2, int threshold);
    int w_sub(string caractere_grafo, string caractere_sequence);
   

public:
    SequenceGraph sequenceGraph;
    my_Bifrost(int k);
    my_Bifrost(int k, string nomeArquivo);
    CompactedDBG<myData> newDeBruijnGraph(int k, string nomeArquivo, bool detalhes);
    void clearMarking(CompactedDBG<myData>& ccdbg);
    void insertSequence(string sequence);
    void clear();
    int size();
    bool haskmer(const Kmer kmer);
    vector<int> findAnchors(string sequence);
    bool compareKmersWithGraph(string kmer, int threshold);
    void BFS_Iterative(const UnitigMap<myData>& ucm);
    void BFS_Iterative(const UnitigMap<myData>& ucm, int distance);
    void BFS_IterativeByDistanceAndInsertKmers(const UnitigMap<myData>& ucm,  my_Bifrost &sub_graph, int distance, bool sucessor);
    bool BFS_Iterative(const UnitigMap<myData>& ucm, string kmer, int threshold);
    bool BFS_IterativeAndInsertKmers(const UnitigMap<myData>& ucm, my_Bifrost &sub_graph, string kmer, int threshold, bool sucessor);
    void traverse();
    void traverse(int distance);
    bool traverse(string kmer, int threshold);
    void traverseByDistance(my_Bifrost &sub_graph, string kmer, int distance, bool sucessor);
    void insertSpecialsKmers();
    void convertDbgToSequenceGraph();
    pair<SequenceGraph, int>  buildMultilayerGraph(SequenceGraph grafo, string sequence);
    void shortestPath(SequenceGraph grafo, int src, int dest, int W);
    pair<vector<pair<int,string>>, int> dijkstra(SequenceGraph grafo, int orig, int dest);
    string verifyEdge(int u, int v, int tamGraph);
    pair<list<string>, string> buildMapping(vector<pair<int,string>> retorno, my_Bifrost dbg_gap, SequenceGraph graph);
    pair<int, int> vertice_inicial_final(string cabeca, string cauda);
    string findKmerByIndex(int index);

};

my_Bifrost::my_Bifrost(int k)
{
    this->k = k;
    this->cdbg = CompactedDBG<myData>(k);
}

my_Bifrost::my_Bifrost(int k, string nomeArquivo)
{
    this->k = k;
    this->cdbg = this->newDeBruijnGraph(k, nomeArquivo, false);
}

CompactedDBG<myData> my_Bifrost::newDeBruijnGraph(int k, string nomeArquivo, bool detalhes = false)
{
    string linha;
	fstream meuArquivo;
	meuArquivo.open(nomeArquivo);
    CompactedDBG<myData> cdbg(k);
    bool status;

    //cdbg.read(nomeArquivo);
    //cout << "Qtd. de kmers no grafo " << cdbg.nbKmers() << endl;
    //return cdbg;


	if (!meuArquivo) {
		cout << "Arquivo " << nomeArquivo << " de kmers não encontrado" << endl;
	}
	else {
        if (detalhes) cout << "Criando o Grafo" << endl;
        getline(meuArquivo, linha);
        while (getline(meuArquivo, linha))
        {
            if (detalhes) cout << "Adicionando " << linha << endl;
            cdbg.add(linha, detalhes);
            getline(meuArquivo, linha);
            if (detalhes) cout << "trash " << linha << endl;
        }
		meuArquivo.close();
    }    
    return cdbg;
}

void my_Bifrost::insertSequence(string sequence)
{
    // cout << "insere " << sequence << endl;
    this->cdbg.add(sequence);
}

void my_Bifrost::clear()
{
    this->cdbg.clear();
    this->cdbg = CompactedDBG<myData>(this->k);
}

int my_Bifrost::size()
{
    return this->cdbg.nbKmers();
}


bool my_Bifrost::haskmer(const Kmer kmer)
{
    if (this->cdbg.find(kmer).isEmpty)
        return 0;    
    return 1;
}

vector<int> my_Bifrost::findAnchors(string sequence)
{
    vector<int> positions;
    string kmer_sequence;
    for (int i = 0; i < sequence.length() - this->k; i++)
    {   
        kmer_sequence = sequence.substr(i, this->k);
        const Kmer km = Kmer(kmer_sequence.c_str()); 
        if (this->haskmer(km))
            positions.push_back(i);
        //else
          //  positions.push_back(0);
    }
    return positions;
}

bool my_Bifrost::compareKmersWithGraph(string kmer, int threshold)
{
    int d = 0;
    for (const auto& unitig : cdbg){
        auto kmer_mapped = unitig.getMappedKmer(0);
        d = this->hammingDistance(kmer_mapped.toString(), kmer, threshold);
        if (d > threshold)
            return 0;   
    }
    return 1;
}

int my_Bifrost::hammingDistance(string kmer1, string kmer2, int threshold)
{
    int i = 0, errors = 0;
    while (i < kmer1.length() && errors <= threshold)
    {
        if (kmer1[i] != kmer2[i])
        {
            errors++;
        }
        i++;
    }
    return errors;
}

void my_Bifrost::BFS_Iterative(const UnitigMap<myData>& ucm)
{
    queue<UnitigMap<myData>> q; // Create queue of unitig to traverse

    UnitigMap<myData> ucm_tmp(ucm); // Create a non-const local copy of unitig given in parameter

    myData* data = ucm_tmp.getData(); // Get DataAccessor from unitig

    if (!data)
    {
        std::cout << "ptr is a null pointer." << endl;
        return;
    }

    data->set_visited(); // Set boolean to indicate unitig was visited
    q.push(ucm_tmp); // Push unitig to traverse on the stack
    
    while (!q.empty()){ // While they are unitigs to traverse in the stack

        ucm_tmp = q.front(); // Get unitig at the front of the queue

        q.pop(); // Delete unitig at the front of the queue

        for (auto& successor : ucm_tmp.getSuccessors()){ // Traverse successors

            myData* data_succ = successor.getData(); // Get boolean from DataAccessor

            if (data_succ->is_not_visited()){ // If boolean indicates the successor was not visited

                data_succ->set_visited(); // Set boolean to indicate successor was visited

                q.push(successor); // Traverse neighbors of successor
            }
        }

        // Traverse predecessors
        for (auto& predecessor : ucm_tmp.getPredecessors()){

            myData* data_pred = predecessor.getData(); // Get DataAccessor from predecessor

            if (data_pred->is_not_visited()){ // If boolean indicates the predecessor was not visited

                data_pred->set_visited(); // Set boolean to indicate predecessor was visited

                q.push(predecessor); // Traverse neighbors of predecessor
            }
        }
    }
}

void my_Bifrost::BFS_Iterative(const UnitigMap<myData>& ucm, int distance)
{
    int distance_tmp = 0;
    queue<UnitigMap<myData>> q; // Create queue of unitig to traverse

    UnitigMap<myData> ucm_tmp(ucm); // Create a non-const local copy of unitig given in parameter

    myData* data = ucm_tmp.getData(); // Get DataAccessor from unitig

    if (!data)
    {
        std::cout << "ptr is a null pointer." << endl;
        return;
    }

    data->set_visited(); // Set boolean to indicate unitig was visited
    q.push(ucm_tmp); // Push unitig to traverse on the stack
    
    while (!q.empty()){ // While they are unitigs to traverse in the stack

        ucm_tmp = q.front(); // Get unitig at the front of the queue

        q.pop(); // Delete unitig at the front of the queue

        distance_tmp = 0;
        for (auto& successor : ucm_tmp.getSuccessors()){ // Traverse successors
            myData* data_succ = successor.getData(); // Get boolean from DataAccessor

            if (data_succ->is_not_visited()){ // If boolean indicates the successor was not visited

                data_succ->set_visited(); // Set boolean to indicate successor was visited

                q.push(successor); // Traverse neighbors of successor
            }
            distance_tmp++;
            if (distance_tmp > distance)
                break;
        }

        // Traverse predecessors
        distance_tmp = 0;
        for (auto& predecessor : ucm_tmp.getPredecessors()){

            myData* data_pred = predecessor.getData(); // Get DataAccessor from predecessor

            if (data_pred->is_not_visited()){ // If boolean indicates the predecessor was not visited

                data_pred->set_visited(); // Set boolean to indicate predecessor was visited

                q.push(predecessor); // Traverse neighbors of predecessor
            }
            distance_tmp++;
            if (distance_tmp > distance)
                break;
        }
    }
}


void my_Bifrost::BFS_IterativeByDistanceAndInsertKmers(const UnitigMap<myData>& ucm,  my_Bifrost &sub_graph, int distance, bool sucessor)
{
    int distance_tmp = 0;
    queue<UnitigMap<myData>> q; // Create queue of unitig to traverse

    UnitigMap<myData> ucm_tmp(ucm); // Create a non-const local copy of unitig given in parameter

    myData* data = ucm_tmp.getData(); // Get DataAccessor from unitig

    if (!data)
    {
        std::cout << "ptr is a null pointer." << endl;
        return;
    }

    data->set_visited(); // Set boolean to indicate unitig was visited
    q.push(ucm_tmp); // Push unitig to traverse on the stack
    distance_tmp = 0;
   
    while (!q.empty()){ // While they are unitigs to traverse in the stack

        if (distance_tmp > distance)
            break;  
        ucm_tmp = q.front(); // Get unitig at the front of the queue

        q.pop(); // Delete unitig at the front of the queue
        if (sucessor == 1)
        {
            // cout << "distance " << distance_tmp << " < " << distance << endl;

            if (ucm_tmp.getSuccessors().hasSuccessors()){
                distance_tmp++;

                for (auto& successor : ucm_tmp.getSuccessors()){ // Traverse successors
                    myData* data_succ = successor.getData(); // Get boolean from DataAccessor

                    if (data_succ->is_not_visited()){ // If boolean indicates the successor was not visited

                        data_succ->set_visited(); // Set boolean to indicate successor was visited

                        q.push(successor); // Traverse neighbors of successor
                        // cout << "Insere ---> " << successor.getMappedKmer(0).toString() << endl;
                        sub_graph.insertSequence(successor.getMappedKmer(0).toString()); // colocando kmers no subgraph
                    }
                    distance_tmp++;
                }
            }
        }
        else
        { // Traverse predecessors 
            // cout << "distance " << distance_tmp << " < " << distance << endl;

            if (ucm_tmp.getPredecessors().hasPredecessors())
            {
                distance_tmp++;
                for (auto& predecessor : ucm_tmp.getPredecessors()){

                    myData* data_pred = predecessor.getData(); // Get DataAccessor from predecessor
   
                    if (data_pred->is_not_visited()){ // If boolean indicates the predecessor was not visited

                        data_pred->set_visited(); // Set boolean to indicate predecessor was visited

                        q.push(predecessor); // Traverse neighbors of predecessor
                        // cout << "Insere <---- " << predecessor.getMappedKmer(0).toString() << endl;
                        sub_graph.insertSequence(predecessor.getMappedKmer(0).toString()); // colocando kmers no subgraph
                    }

                }
            }
        }
    }
}

bool my_Bifrost::BFS_IterativeAndInsertKmers(const UnitigMap<myData>& ucm, my_Bifrost &sub_graph, string kmer, int threshold, bool sucessor)
{
    queue<UnitigMap<myData>> q; // Create queue of unitig to traverse

    UnitigMap<myData> ucm_tmp(ucm); // Create a non-const local copy of unitig given in parameter

    myData* data = ucm_tmp.getData(); // Get DataAccessor from unitig

    if (!data)
    {
        std::cout << "ptr is a null pointer." << endl;
        return 0;
    }

    data->set_visited(); // Set boolean to indicate unitig was visited
    q.push(ucm_tmp); // Push unitig to traverse on the stack
    
    while (!q.empty()){ // While they are unitigs to traverse in the stack

        ucm_tmp = q.front(); // Get unitig at the front of the queue

        int d = this->hammingDistance(ucm_tmp.getMappedKmer(0).toString(), kmer, threshold);
        if (d > threshold)
            return 0;

        q.pop(); // Delete unitig at the front of the queue

        if (sucessor == 1)
        {
            cout << ucm_tmp.getMappedKmer(0).toString() << " tem sucessor? " << ucm_tmp.getSuccessors().hasSuccessors() << endl;
            for (auto& successor : ucm_tmp.getSuccessors()){ // Traverse successors

                myData* data_succ = successor.getData(); // Get boolean from DataAccessor

                if (data_succ->is_not_visited()){ // If boolean indicates the successor was not visited

                    data_succ->set_visited(); // Set boolean to indicate successor was visited

                    q.push(successor); // Traverse neighbors of successor

                    cout << "insere ---> " << successor.getMappedKmer(0).toString() << endl;
                    sub_graph.insertSequence(successor.getMappedKmer(0).toString()); // colocando kmers no subgraph
                    
                }
            }
        }else {
            // Traverse predecessors
            for (auto& predecessor : ucm_tmp.getPredecessors()){

                myData* data_pred = predecessor.getData(); // Get DataAccessor from predecessor

                if (data_pred->is_not_visited()){ // If boolean indicates the predecessor was not visited

                    data_pred->set_visited(); // Set boolean to indicate predecessor was visited

                    q.push(predecessor); // Traverse neighbors of predecessor
        
                    cout << "insere ---> " << predecessor.getMappedKmer(0).toString() << endl;
                    sub_graph.insertSequence(predecessor.getMappedKmer(0).toString()); // colocando kmers no subgraph
                    
                }
            }
        }
    }
    return 1;
}

bool my_Bifrost::BFS_Iterative(const UnitigMap<myData>& ucm, string kmer, int threshold)
{
    queue<UnitigMap<myData>> q; // Create queue of unitig to traverse

    UnitigMap<myData> ucm_tmp(ucm); // Create a non-const local copy of unitig given in parameter

    myData* data = ucm_tmp.getData(); // Get DataAccessor from unitig

    if (!data)
    {
        std::cout << "ptr is a null pointer." << endl;
        return 0;
    }

    data->set_visited(); // Set boolean to indicate unitig was visited
    q.push(ucm_tmp); // Push unitig to traverse on the stack
    
    while (!q.empty()){ // While they are unitigs to traverse in the stack

        ucm_tmp = q.front(); // Get unitig at the front of the queue

        int d = this->hammingDistance(ucm_tmp.getMappedKmer(0).toString(), kmer, threshold);
        if (d > threshold)
            return 0;

        q.pop(); // Delete unitig at the front of the queue

        for (auto& successor : ucm_tmp.getSuccessors()){ // Traverse successors

            myData* data_succ = successor.getData(); // Get boolean from DataAccessor

            if (data_succ->is_not_visited()){ // If boolean indicates the successor was not visited

                data_succ->set_visited(); // Set boolean to indicate successor was visited

                q.push(successor); // Traverse neighbors of successor
            }
        }

        // Traverse predecessors
        for (auto& predecessor : ucm_tmp.getPredecessors()){

            myData* data_pred = predecessor.getData(); // Get DataAccessor from predecessor

            if (data_pred->is_not_visited()){ // If boolean indicates the predecessor was not visited

                data_pred->set_visited(); // Set boolean to indicate predecessor was visited

                q.push(predecessor); // Traverse neighbors of predecessor
            }
        }
    }

    return 1;
}

void my_Bifrost::traverse()
{
    for (auto& unitig : this->cdbg){ // Iterate over unitigs of a colored de Bruijn graph
        BFS_Iterative(unitig); // Traverse neighbors of unitig in an iterative manner
    }
}

void my_Bifrost::traverse(int distance)
{
    for (auto& unitig : this->cdbg){ // Iterate over unitigs of a colored de Bruijn graph
        BFS_Iterative(unitig, distance); // Traverse neighbors of unitig in an iterative manner
    }
}

bool my_Bifrost::traverse(string kmer, int threshold)
{
    for (auto& unitig : this->cdbg){ // Iterate over unitigs of a colored de Bruijn graph
        if (!BFS_Iterative(unitig, kmer, threshold)) // Traverse neighbors of unitig in an iterative manner
            return 0;
    }
    return 1;
}

void my_Bifrost::traverseByDistance(my_Bifrost &sub_graph, string kmer, int distance, bool sucessor)
{
    const Kmer km = Kmer(kmer.c_str());
    auto unitig = this->cdbg.find(km);
    BFS_IterativeByDistanceAndInsertKmers(unitig, sub_graph, distance, sucessor);    
    // cout << "atravessei o gap " << distance << endl;
}

void my_Bifrost::clearMarking(CompactedDBG<myData>& ccdbg){

    for (const auto& unitig : ccdbg)
    {
        myData* data_ucm = unitig.getData();
        data_ucm->clear(unitig);
    }
}

void my_Bifrost::insertSpecialsKmers()
{
    unordered_map<string, list<string>>:: iterator itr;
    string bases[] = {"A", "C", "G", "T"}, kmer, kmer_aux; int val;
    list<string>::iterator it; list<string> bkp;

    for (auto unitig : this->cdbg)
    {
        if (!unitig.getPredecessors().hasPredecessors())
        {
            string kmer = unitig.referenceUnitigToString().substr(0, unitig.dist + k - 1); // prefixo com k-1
            for (int i = 1; i < kmer.length(); i++)
            {
                kmer_aux = "";
                for (int j = 0; j < i; j++)
                    kmer_aux = kmer_aux + "$";                
                kmer_aux = kmer_aux + kmer.substr(0,this->k-i);
                //kmerSpecialAndKmer[kmer_aux] = kmer;
                this->cdbg.add(kmer_aux);
            }            
        }
    }      
}

// pair<unordered_map<int, string>, SequenceGraph>
void my_Bifrost::convertDbgToSequenceGraph()
{
    unordered_map<int, string> kmerAndNode;
    int qtdNode = 0;
    for (auto& unitig : this->cdbg)
    {
        auto unitig_aux = unitig.mappedSequenceToString();
        myData* data = unitig.getData(); // Get DataAccessor from unitig
        data->set_initial_node(qtdNode);
        data->set_node_in_sequence_graph(qtdNode);
        qtdNode = qtdNode + unitig_aux.length();
    }

    sequenceGraph.initilizeSequenceGraph(qtdNode, this->k);

    for (auto& unitig : this->cdbg)
    {
        auto unitig_aux = unitig.mappedSequenceToString();
        myData* data = unitig.getData(); // Get DataAccessor from unitig
        int initial = data->get_initial_node();

        string kmer = unitig_aux.substr(0, k);
        for (int i = 0; i < unitig_aux.length(); i++)
        {
            if (i % this->k == 0 && i < unitig_aux.length() - (this->k - 1))
                kmer = unitig_aux.substr(i, k);
            if(i == 0)
                sequenceGraph.alterarValorVerticeInicial(initial + i, 1);
            sequenceGraph.insertNode(initial + i, unitig_aux.substr(i, 1));
            // kmerAndNode[initial + i] = kmer; // mapeando kmer e node
        }
        for (int i = 0; i < unitig_aux.length() - 1; i++)
        {
            sequenceGraph.insertEdge(initial + i, initial + i + 1, 0);
        }
    } 

    for (auto& unitig : this->cdbg)
    {
        auto unitig_aux = unitig.mappedSequenceToString();
        myData* data = unitig.getData(); // Get DataAccessor from unitig
        int initial = data->get_initial_node();

        if (unitig.getSuccessors().hasSuccessors())
        {
            int nodeUnitig = initial + unitig_aux.length() - 1;
            for (auto sucessor : unitig.getSuccessors())
            {
                myData* data = sucessor.getData(); // Get DataAccessor from unitig
                int nodeSucessor = data->get_initial_node() + sucessor.mappedSequenceToString().length() - 1;
                sequenceGraph.insertEdge(nodeUnitig, nodeSucessor, 0);
            }
        }
    }
}

pair<SequenceGraph, int> my_Bifrost::buildMultilayerGraph(SequenceGraph grafo, string sequence)
{
    int V = grafo.getV(), m = sequence.length(), vertice_atual = 0, vertice_inicial = 0, vertice_final = -1, vertice_atual_aux, controle;
    int m_v = m * (V + 1) + 2; // quantidade de vertice do grafo multicamadas
    SequenceGraph m_grafo(m_v, grafo.getK());
    int *mapeamento;

    mapeamento = new (nothrow) int[V + 1];
    sequenceGraphAndMulticamada = new (nothrow) vector<int>[m_v];  
    if (mapeamento == nullptr || sequenceGraphAndMulticamada == nullptr)
    {
        cerr << "error allocation multlayer graph" << endl;
        return make_pair(m_grafo, vertice_final);
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
                    {
                        m_grafo.insertEdge(vertice_atual - (V + 1), mapeamento[j], w_sub(grafo.getBase(j), sequence.substr(i-1,1)));                  
                    }
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
    //this->insertInitialAndEndNode(vertice_inicial, vertice_final);

    // deletando vetores sem utilizacao
    delete [] mapeamento;
    return make_pair(m_grafo, vertice_final);
}

int my_Bifrost::w_sub(string caractere_grafo, string caractere_sequence)
{
    if (caractere_grafo.compare(caractere_sequence) == 0)
        return 0;
    return sub;
}

// Prints shortest paths from src to all other vertices.
// W is the maximum weight of an edge
void my_Bifrost::shortestPath(SequenceGraph grafo, int src, int dest, int W)
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

// Dijkstra
pair<vector<pair<int,string>>, int> my_Bifrost::dijkstra(SequenceGraph grafo, int orig, int dest)
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

string my_Bifrost::verifyEdge(int u, int v, int tamGraph)
{
    int lim = u + (tamGraph + 1);

    if (v == lim)
        return "del";
    if (v >= lim - ((tamGraph+1)/2) && v < lim)
        return "sub";
    if (v > lim)
        return "sub";
    return "ins";
}

pair<list<string>, string> my_Bifrost::buildMapping(vector<pair<int,string>> retorno, my_Bifrost dbg_gap, SequenceGraph graph)
{
    int primeiro = 0, indice, anterior = 0, details = 0, k = graph.getK(), kmer_count = 0;
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

            tmp = this->verifyEdge((*it).first, anterior, graph.getV());
            if (details == 1)
                cout << "(" << tmp << ") ";
            anterior = (*it).first;    
            indice = this->sequenceGraphAndMulticamada[(*it).first].front(); 

            if (indice != -1)
            {   
                auto kmer = dbg_gap.findKmerByIndex(indice);      
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
    return  make_pair(kmers,aux.substr(0, aux.length() - 1));
}

pair<int, int> my_Bifrost::vertice_inicial_final(string cabeca, string cauda)
{
	int c = 0, t = 0; 
    const Kmer km_cabeca = Kmer(cabeca.c_str());
    const Kmer km_cauda = Kmer(cauda.c_str());

    if(!cdbg.find(km_cabeca).isEmpty)
    {
        auto unitig = cdbg.find(km_cabeca);
        myData* data_succ = unitig.getData(); 
        c = data_succ->get_node_in_sequence_graph();
    }

    if(!cdbg.find(km_cauda).isEmpty)
    {
        auto unitig = cdbg.find(km_cauda);
        myData* data_succ = unitig.getData(); 
        t = data_succ->get_node_in_sequence_graph();    }

	return make_pair(c,t);
}

string my_Bifrost::findKmerByIndex(int index)
{
    for (const auto& unitig : cdbg){
        myData* data_succ = unitig.getData(); 
        int dif = abs(data_succ->get_node_in_sequence_graph() - index);
        if (dif <= this->k)
        {
            //cout << "kmer teste " << itr->first << endl;
            auto kmer_mapped = unitig.getMappedKmer(0);
            return kmer_mapped.toString();
        }
    }
    return "";
}

/*int main2(int argc, char **argv){

    string nomeArquivo = "/home/lucas/Documentos/sequencias/lordec/reads/one_short_reads.fasta", linha;

    size_t k;
    cin >> k;
    auto cdbg = newDeBruijnGraph(k, nomeArquivo);
    
    cout << "K-mer size is " << cdbg.getK() << endl;
   
    const string sequence = "AACCGAACAGTATCGTGCCATCTTGTATGCCGCGCTCCTG";

    auto retorno = findAnchors(cdbg, sequence);

    for (auto it : retorno)
    {
        cout << it << " ";
    }
    cout << endl;
    return 0;
}*/