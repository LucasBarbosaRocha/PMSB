#include <fstream>
#include <iostream>
#include "../utils/marschallTwoLayers.cpp"
//#include "../utils/myBifrost/my_bifrost.cpp"
#include "../utils/myUtils.cpp"
#include <queue>
#include <bits/stdc++.h>

#define x 0.5
MyUtils utils;
MarschallTwoLayers m2;

int returnLengthGap(int k, int gap)
{
    if (gap > k)
        return gap * x;
    return k;
}

pair<int, int> headMapping(my_Bifrost dbg_gap, int firstPosition, my_Bifrost h, string sequence, int k)
{
    int length, aux, caminho_encontrado = 1;
    list<string> kmer_lista_aux;
    list<int> nodes;
    string resposta;
    int cost = 0;


    if (firstPosition != 0)
    {
        // cout << "HEAD " << endl;
        length = firstPosition * x, aux = 0;

        if (firstPosition > 10)
            length = returnLengthGap(k, firstPosition); 
        else
            length = firstPosition - 0;       
        aux = 0;        
        string kmer_cauda = sequence.substr(firstPosition, k);

        dbg_gap.insertSequence(kmer_cauda);
        h.traverseByDistance(dbg_gap, kmer_cauda, length, 0);

        dbg_gap.convertDbgToSequenceGraph();
        string sequence_aux = sequence.substr(0, k + firstPosition);
        cost = m2.toCalculateTheCost(dbg_gap.sequenceGraph, sequence_aux);

        if (cost == INT_MAX)  
        {
            caminho_encontrado = 0;
            return make_pair(caminho_encontrado, cost);
        }
        // liberando memória
        dbg_gap.sequenceGraph.deleteGraph();
        dbg_gap.clear();

    }
    return make_pair(caminho_encontrado, cost);
}

pair<int, int> internalMapping(my_Bifrost dbg_gap, vector<int> positions, my_Bifrost h, string sequence, int k)
{
    int caminho_encontrado = 0, length, aux, posicao = positions[0];
    string resposta = "", kmer_cabeca = "", kmer_cauda = "", sequence_aux = "";
    list<int> nodes;
    int cost = INT_MAX, cost_aux = 0;
    // cout << "INTERNAL " << endl;
    for(int i = 0; i < positions.size() - 1; i++)
    {     
        caminho_encontrado = 1;

        int dif = positions[i + 1] - positions[i]; 
        if (dif > 1)
        {
            length = returnLengthGap(k, dif); aux = 0;
            kmer_cabeca = sequence.substr(positions[i], k);
            kmer_cauda = sequence.substr(positions[i+1], k);

            // buscando os kmers do gap
            dbg_gap.insertSequence(kmer_cabeca);
            dbg_gap.insertSequence(kmer_cauda);
            
            h.traverseByDistance(dbg_gap, kmer_cabeca, length, 1);
            h.traverseByDistance(dbg_gap, kmer_cauda, length, 0);

            dbg_gap.convertDbgToSequenceGraph();
            sequence_aux = sequence.substr(positions[i], k + dif);
            cost = m2.toCalculateTheCost(dbg_gap.sequenceGraph, sequence_aux);

            if (cost == INT_MAX)
            {
                i = positions.size() + 1;
                caminho_encontrado = positions[i + 1];
                break;
            } else
            {
                //cout << "achei " << cost << endl;
                cost_aux += cost;
            }
            // limpando memória
            dbg_gap.sequenceGraph.deleteGraph();
            dbg_gap.clear();
        }else
            posicao++;
    }
    return make_pair(caminho_encontrado, cost_aux);
}

pair<int, int> tailMapping(my_Bifrost dbg_gap, int lastPosition, my_Bifrost h, string sequence, int k)
{
    int length = (sequence.length() - lastPosition) * x, aux = 0, caminho_encontrado = 0, dif;
    int cost = 0;
    list<int> nodes;
    string resposta = "", kmer_cabeca = "", sequence_aux = "";

    if (sequence.length() - lastPosition >= k)
    {
        length = returnLengthGap(k, sequence.length() - lastPosition);
        kmer_cabeca = sequence.substr(lastPosition, k);

        dbg_gap.insertSequence(kmer_cabeca);
        h.traverseByDistance(dbg_gap, kmer_cabeca, length, 1);

        dbg_gap.convertDbgToSequenceGraph();
        dif = sequence.length() - lastPosition;    
        sequence_aux = sequence.substr(lastPosition, k + dif);
        cost = m2.toCalculateTheCost(dbg_gap.sequenceGraph, sequence_aux);

        if (cost != INT_MAX)      
        {
            caminho_encontrado = 1; 
        }     
        dbg_gap.sequenceGraph.deleteGraph();
        dbg_gap.clear();  
    }
    return make_pair(caminho_encontrado, cost);
}

int mapeamento(my_Bifrost h, string sequence, int k)
{
    int posicao = 0, length, aux, caminho_encontrado = 1;
    string resposta = "";
    vector<int> posicoesValidas;
    list<string> kmer_lista_aux;
    my_Bifrost dbg_gap(k);
    pair<int, int> status;
    int cost_aux = sequence.length();


    posicoesValidas = h.findAnchors(sequence);


    cout << "Qtd. Anchros " << posicoesValidas.size() << endl;

    if (posicoesValidas.size() > 0)
    {
        status = headMapping(dbg_gap, posicoesValidas[0], h, sequence, k);
        
        if (status.second >= INT_MAX)
        {   
            int extra_cost = sequence.length() - posicoesValidas[0];
            return status.second + extra_cost;
        }

        // cout << "Head path " << status.first << endl;

        cost_aux = status.second;

        status = internalMapping(dbg_gap, posicoesValidas, h, sequence, k);

        if (status.second >= INT_MAX)
        {
            int extra_cost = sequence.length() - status.first;
            return cost_aux + status.second + extra_cost;
        }

        // cout << "Internal path " << status.first << endl;

        cost_aux += status.second;

        status = tailMapping(dbg_gap, posicoesValidas[posicoesValidas.size() - 1], h, sequence, k);

        if (status.second >= INT_MAX)
        {
            int extra_cost = sequence.length() - posicoesValidas[posicoesValidas.size() - 1];
            return cost_aux + status.second + extra_cost;
        }
        cost_aux += status.second;
       
        // cout << "Tail path " << status.first << endl;
    } 
    return cost_aux;
}

int main(int argc, char *argv[])
{
    string line;
    if(utils.verifyData(argc, argv) == 1)
        exit (0); 
        
    ifstream file(utils.nameSequenceArchive);
    my_Bifrost bf(utils.k, utils.nameArchive);  

    while(getline(file, line))
    {
        //utils.readSequence(utils.nameSequenceArchive);  
        getline(file, line);
        cout << "Size L.Read " << line.size() << endl;
        utils.sequence = line;
        // mapeamento
        auto retorno = mapeamento(bf, utils.sequence, utils.k); 
        cout << retorno << endl << endl;
    }     
    return 0;
}
