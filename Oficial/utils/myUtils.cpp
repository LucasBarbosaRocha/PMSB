#include <string>
#include <iostream>

using namespace std;

class MyUtils
{
public:
    string sequence;
    string nameSequenceArchive;
    string nameArchive;
    int k;
    int typeGraph;

    MyUtils(){};
    int verifyData(int argc, char *argv[]);
    void readSequence(string archiveName);
};

int MyUtils::verifyData(int argc, char *argv[])
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
            cout << "-s nameSequenceArchive -g nameGraphArchive -k kmerSize -t typeSequenceGraph (0 - traditional, 1 - simplified) " << endl;
        return 1;
    }

    if (argc == 9)
    {
        aux = argv[1];
        if (aux.compare("-s") == 0)
            nameSequenceArchive = argv[2];
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
        aux = argv[7];
        if (aux.compare("-t") == 0)
            typeGraph = atoi(argv[8]);
        else
        {
            cout << "Error: digite -help" << endl;
            return 1;
        }
        return 0;
    }

    cout << "Error: digite -help" << endl;
    return 1;
}

void MyUtils::readSequence(string archiveName)
{
    string line;
	fstream myArchive;
	myArchive.open(archiveName);

	if (!myArchive) {
		cout << "Arquivo " << archiveName << " de sequence não encontrado" << endl;
	}
	else {
        getline(myArchive, line);     
        getline(myArchive, line);     
		myArchive.close();
    }    
    this->sequence = line;
}