// de Bruijn sequence mapping Tool with graph Label change
#include <iostream>
#include "Hungarian.h"
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"
#include <string>

MyUtils utils;

using namespace std;

vector< vector<double> > createCostMatrix(Hash h, string sequence, int k)
{
	string kmer; vector< vector<double>> costMatrix;
	for (int i = 0; i <= sequence.length() - k; i++)
	{
		kmer = sequence.substr(i, k);
		auto costs = h.getCosts(kmer);
		costMatrix.push_back(costs);
	}

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
	unordered_map<string, int> kmersInSequence;
	string kmer;

    h.populateGraph(utils.nameArchive, false);    
	h.enrichGraph();   

    while(getline(file, line))
    {
        //utils.readSequence(utils.nameSequenceArchive);  
        getline(file, line);
		utils.sequence = line;
        cout << "Kmers in Gk: " << h.getQtdKmers() << endl;
        cout << "Size L.Read " << utils.sequence.size() << endl;
		for (int i = 0; i < utils.sequence.size() - utils.k; i++)
		{
			kmer = utils.sequence.substr(i, utils.k);
			kmersInSequence[kmer] = 1;
		}
		cout << "qtd kmers in s " << kmersInSequence.size() << endl;

		if (h.getQtdKmers() >= kmersInSequence.size())
		{
			kmersInSequence.clear();	
			// mapeamento
			auto matrix = createCostMatrix(h, utils.sequence, utils.k);
			/*for (int i = 0; i < matrix.size(); i++) {
					for (int j = 0; j < matrix[i].size(); j++)
						cout << matrix[i][j] << " ";
					cout << endl;
				}*/
			double cost = HungAlgo.Solve(matrix, assignment);
			for (unsigned int x = 0; x < matrix.size(); x++)
			{
				maping[assignment[x]] = x;
			}
			
			cout << "cost: " << cost << endl;

			h.changeGraph(maping, utils.sequence);

			int count = 0;
			for (int i = 0; i < utils.sequence.size() - ((utils.k) - 1); i++)
			{
				string kmer = utils.sequence.substr(i, utils.k);
				if(h.contains(kmer))
				{
					count++;
				}
			}

			cout << "kmers after maping: " << count << endl;
		} else {
			cout << "Qty. kmers (" << h.getQtdKmers() << ") in Gk need to be >= kmers in sequence to pair" << endl;
		}
	}     
    return 0;
}
