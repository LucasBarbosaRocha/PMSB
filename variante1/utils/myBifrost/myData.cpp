#include <string>
#include <stack>
#include <queue>
#include <unordered_set>

#include "bifrost/src/CompactedDBG.hpp"

using namespace std;

class myData  : CDBG_Data_t<myData>
{
    private:
        const static uint8_t NOT_VISITED_SEEN = 0x0;
        const static uint8_t VISITED = 0x1;
        const static uint8_t SEEN = 0x2;
        uint8_t b;
        int initial_node;
        int node_in_sequence_graph;
        unordered_map<int, string> map_node_kmer;

    public:
        // Clear method for CompactedDBG
        void clear(const UnitigMap<myData>& um_dest){

        	set_not_seen_visited(); // Set the new unitig to "not seen nor visited"
        }

        // Concatenation method for ColoredCDBG
        void concat(const UnitigMap<myData>& um_dest, const UnitigMap<myData>& um_src){

            // When concatenating the reference unitig of um_src to the reference unitig of um_dest,
            // we set the boolean of the new unitig to "not visited" because it will be a new unitig.

            set_not_seen_visited(); // Set the new unitig to "not seen nor visited"
        }

        // Extraction method for ColoredCDBG
        void extract(const UnitigMap<myData>& um_src, bool last_extraction) {

            // This function creates a new unitig which is a sub-unitig from the reference unitig of um_src.
            // The new unitig created is set to "not seen nor visited" as a measure of precaution (it is already
            // initiated by default to "not visited" in the constructor)

            set_not_seen_visited(); // Set the new unitig to "not seen nor visited"
        }

        // Methods myData::merge() and myData

        string toString() const {

            if (is_visited()) return string("visited");
            if (is_seen()) return string("seen");

            return string("Not seen nor visited");
        }

        inline void set_visited() { b = VISITED; } // Set the boolean to "visited"
        inline void set_seen() { b = SEEN; } // Set the boolean to "seen"
        inline void set_not_seen_visited() { b = NOT_VISITED_SEEN; } // Set the boolean to "not seen and not visited"
        inline void set_initial_node(int node) { this->initial_node = node; }
        inline int get_initial_node() { return this->initial_node; }

        inline void set_node_in_sequence_graph(int node_in_sequence_graph) { this->node_in_sequence_graph = node_in_sequence_graph; }
        inline int get_node_in_sequence_graph() { return this->node_in_sequence_graph; }

        inline bool is_visited() const { return (b == VISITED); } // return if the boolean is "visited"
        inline bool is_not_visited() const { return !is_visited(); } // return if the boolean is "not visited"

        inline bool is_seen() const { return (b == SEEN); } // return if the boolean is "seen"
        inline bool is_not_seen() const { return !is_seen(); } // return if the boolean is "not seen"

        string getKmerBySequenceGraphNode(int node)
        { return this->map_node_kmer[node]; }

};

/*

void BFS_Iterative(const UnitigMap<myData>& ucm)
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

        cout << ucm_tmp.getMappedKmer(0).toString() << " ";

        for (auto& successor : ucm_tmp.getSuccessors()){ // Traverse successors

            myData* data_succ = successor.getData(); // Get boolean from DataAccessor
            cout << successor.getMappedKmer(0).toString() << " ";

            if (data_succ->is_not_visited()){ // If boolean indicates the successor was not visited

                data_succ->set_visited(); // Set boolean to indicate successor was visited

                q.push(successor); // Traverse neighbors of successor
            }
        }

        cout << endl;

        // Traverse predecessors
        for (auto& predecessor : ucm_tmp.getPredecessors()){

            // Get DataAccessor from predecessor
            myData* data_pred = predecessor.getData(); // Get boolean from DataAccessor

            if (data_pred->is_not_visited()){ // If boolean indicates the predecessor was not visited

                data_pred->set_visited(); // Set boolean to indicate predecessor was visited

                q.push(predecessor); // Traverse neighbors of predecessor
            }
        }
    }
}

void traverse(CompactedDBG<myData> cdbg)
{
    for (auto& unitig : cdbg){ // Iterate over unitigs of a colored de Bruijn graph
        BFS_Iterative(unitig); // Traverse neighbors of unitig in an iterative manner
    }
    //clearMarking(this->cdbg);
}


CompactedDBG<myData> newDeBruijnGraph(int k, string nomeArquivo, bool detalhes = false)
{
    string linha;
	fstream meuArquivo;
	meuArquivo.open(nomeArquivo);
    CompactedDBG<myData> cdbg(k);

    //cdbg.read(nomeArquivo, "cores.txt");

	if (!meuArquivo) {
		cout << "Arquivo " << nomeArquivo << " de kmers não encontrado" << endl;
	}
	else {
        if (detalhes) cout << "Criando o Grafo" << endl;
        getline(meuArquivo, linha);
        while (getline(meuArquivo, linha))
        {
            if (detalhes) cout << "Adicionando " << linha << endl;
            cdbg.add(linha); 
            getline(meuArquivo, linha);
        }
		meuArquivo.close();
        cout << "Grafo criado." << endl;
    }    
    return cdbg;
}

int main(int argc, char **argv)
{    
    string nomeArquivo = "/home/lucas/Documentos/sequencias/lordec/reads/one_long_read.fa", linha;
    size_t k;
    cin >> k;
    int threshold = k * 0.5;
    const string sequence = "AACCGAACAGTATCGTGCCATCTTGTATGCCGCGCTCCTG";
    CompactedDBG<myData> graph = newDeBruijnGraph(k, nomeArquivo); 

    traverse(graph);    
    
    return 0;
}*/