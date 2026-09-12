// de Bruijn sequence mapping Tool with graph Label change
#include <iostream>
#include "Hungarian.h"
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"
#include <string>

MyUtils utils;

using namespace std;

vector< vector<double> > createCostMatrix(Hash &h, const vector<string> &kmers)
{
	vector< vector<double>> costMatrix;
	for (const auto &kmer : kmers)
		costMatrix.push_back(h.getCosts(kmer));

	return costMatrix;
}

int main(int argc, char *argv[])
{
    string line;
    if(utils.verifyData(argc, argv) == 1)
        exit (0); 
        
    Hash h(utils.k);   
    ifstream file(utils.nameSequenceArchive);
	HungarianAlgorithm HungAlgo;
	vector<int> assignment;
	unordered_map<int, int> maping;
	string kmer;

    h.populateGraph(utils.nameArchive, false);
	h.enrichGraph();

    while(getline(file, line))
    {
        //utils.readSequence(utils.nameSequenceArchive);
        getline(file, line);
		utils.sequence = line;
        cout << "Kmers no grafo: " << h.getQtdKmers() << endl;
        cout << "Comprimento da sequência: " << utils.sequence.size() << endl;

		// k(s): conjunto dos k-mers DISTINTOS de s (Cap. 6 da tese), guardando
		// a posição da primeira ocorrência de cada um para reconstruir o
		// k-mer em changeGraph.
		vector<string> distinctKmers;
		vector<int> distinctPositions;
		unordered_map<string, int> seen;
		if (utils.sequence.length() >= (size_t)utils.k)
		{
			for (size_t i = 0; i <= utils.sequence.length() - utils.k; i++)
			{
				kmer = utils.sequence.substr(i, utils.k);
				if (seen.find(kmer) == seen.end())
				{
					seen[kmer] = distinctKmers.size();
					distinctKmers.push_back(kmer);
					distinctPositions.push_back((int)i);
				}
			}
		}
		cout << "Quantidade de Kmers na sequência presente no grafo (antes do mapeamento): " << distinctKmers.size() << endl;

		if (h.getQtdKmers() >= distinctKmers.size())
		{
			// mapeamento
			auto matrix = createCostMatrix(h, distinctKmers);
			/*for (int i = 0; i < matrix.size(); i++) {
					for (int j = 0; j < matrix[i].size(); j++)
						cout << matrix[i][j] << " ";
					cout << endl;
				}*/
			double cost = HungAlgo.Solve(matrix, assignment);
			for (unsigned int x = 0; x < matrix.size(); x++)
			{
				maping[assignment[x]] = distinctPositions[x];
			}

			cout << "Custo: " << cost << endl;

			h.changeGraph(maping, utils.sequence);

			int count = 0;
			if (utils.sequence.length() >= (size_t)utils.k)
			{
				for (size_t i = 0; i <= utils.sequence.length() - utils.k; i++)
				{
					string kmer = utils.sequence.substr(i, utils.k);
					if(h.contains(kmer))
					{
						count++;
					}
				}
			}

			cout << "Quantidade de Kmers na sequência presente no grafo (após o mapeamento): " << count << endl;
		} else {
			cout << "Qty. kmers (" << h.getQtdKmers() << ") in Gk need to be >= kmers in sequence to pair" << endl;
		}
	}
    return 0;
}
