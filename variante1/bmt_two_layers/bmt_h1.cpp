#include <fstream>
#include <iostream>
#include "../utils/marschallTwoLayers.cpp"
#include "../utils/myUtils.cpp"
#include <queue>
#include <bits/stdc++.h>

#define x 0.5
MyUtils utils;
MarschallTwoLayers m2;

pair<int, int> headMapping(int firstPosition, Hash h, string sequence, int k)
{
    Hash dbg_gap(k);
    int length, aux, caminho_encontrado = 1;
    list<string> kmer_lista_aux;
    list<int> nodes;
    string resposta;
    MarschallTwoLayers m2;
    int cost = 0;

    if (firstPosition != 0)
    {
        // cout << "HEAD " << endl;
        length = firstPosition * x, aux = 0;
        string kmer_cauda = sequence.substr(firstPosition, k);

        // buscando os kmers do gap
        int final_aux = h.verifyLabelExistisAndReturnVerticeIndex(kmer_cauda);
    
        if (final_aux == -1)
        {
            caminho_encontrado = 0; resposta = "";
            return make_pair(caminho_encontrado, cost);
        }
        final_aux += (k-1);

        nodes = h.sequenceGraphReverse.bfs(final_aux, length);
        dbg_gap.insertKmer(kmer_cauda);        
        h.insertKmersByNodes(nodes, dbg_gap);

        // TODO CONFERIR os nós do do grafo simplificado
        if (utils.typeGraph == 0)
            dbg_gap.dbgToTraditionalSequenceGraph(0, 0); 
        else
            dbg_gap.dbgToSimplifiedSequenceGraph(0); 

        string sequence_aux = sequence.substr(0, k + firstPosition);
        cost = m2.toCalculateTheCost(dbg_gap.sequenceGraph, sequence_aux);

        if (cost == INT_MAX)  
        {
            caminho_encontrado = 0;
            return make_pair(caminho_encontrado, cost);
        }
        // liberando memória
    }
    return make_pair(caminho_encontrado, cost);
}

pair<int, int> internalMapping(vector<int> positions, Hash h, string sequence, int k)
{
    int caminho_encontrado = 0, length, aux, posicao = positions[0];
    string resposta;
    list<int> nodes;
    Hash dbg_gap(k);
    MarschallTwoLayers m2;
    int cost = INT_MAX, cost_aux = 0;
    // cout << "INTERNAL " << endl;
    for(int i = 0; i < positions.size() - 1; i++)
    {     
        caminho_encontrado = 1;

        int dif = positions[i + 1] - positions[i]; 
        if (dif > 1)
        {
            length = dif * x, aux = 0;
            string kmer_cabeca = sequence.substr(positions[i], k);
            string kmer_cauda = sequence.substr(positions[i+1], k);
            // buscando os kmers do gap
            int inicial_aux = h.verifyLabelExistisAndReturnVerticeIndex(kmer_cabeca);
            int final_aux = h.verifyLabelExistisAndReturnVerticeIndex(kmer_cauda);
        
            if (final_aux == -1 || inicial_aux == -1)
            {
                caminho_encontrado = positions[i + 1]; resposta = "";
                return make_pair(caminho_encontrado, cost);
            }

            final_aux += (k-1);

            dbg_gap.deleteHash();
            dbg_gap.insertKmer(kmer_cabeca);
            dbg_gap.insertKmer(kmer_cauda);
            nodes = h.sequenceGraph.bfs(inicial_aux, length);
            //cout << "1-nodes " << nodes.size() << endl;
            h.insertKmersByNodes(nodes, dbg_gap);
            nodes = h.sequenceGraphReverse.bfs(final_aux, length);
            //cout << "2-nodes " << nodes.size() << endl;
            h.insertKmersByNodes(nodes, dbg_gap);

            //cout << "nodes " << dbg_gap.getQtdKmers() << endl;

            if (utils.typeGraph == 0)
                dbg_gap.dbgToTraditionalSequenceGraph(0, 0);        
            else
                dbg_gap.dbgToSimplifiedSequenceGraph(0);

            string sequence_aux = sequence.substr(positions[i], k + dif);
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
        }else
            posicao++;
        
    }
    return make_pair(caminho_encontrado, cost_aux);
}

pair<int, int> tailMapping(int lastPosition, Hash h, string sequence, int k)
{
    int length = (sequence.length() - lastPosition) * x, aux = 0, caminho_encontrado = 0;
    Hash dbg_gap(k);
    MarschallTwoLayers m2;
    int cost = 0;
    list<int> nodes;
    string resposta = "";

    if (sequence.length() - lastPosition >= k)
    {
        // cout << "TAIL " << endl;

        string kmer_cabeca = sequence.substr(lastPosition, k);
        dbg_gap.insertKmer(kmer_cabeca);   

        // buscando os kmers do gap
        int inicial_aux = h.verifyLabelExistisAndReturnVerticeIndex(kmer_cabeca);
    
        if (inicial_aux == -1)
        {
            caminho_encontrado = 0; resposta = "";
            return make_pair(caminho_encontrado, cost);
        }
        dbg_gap.deleteHash();
        dbg_gap.insertKmer(kmer_cabeca);
        nodes = h.sequenceGraph.bfs(inicial_aux, length);
        h.insertKmersByNodes(nodes, dbg_gap);

        if (utils.typeGraph == 0)
            dbg_gap.dbgToTraditionalSequenceGraph(0, 0);
        else
            dbg_gap.dbgToSimplifiedSequenceGraph(0);

        int dif = sequence.length() - lastPosition;    
        string sequence_aux = sequence.substr(lastPosition, k + dif);
        cost = m2.toCalculateTheCost(dbg_gap.sequenceGraph, sequence_aux);

        if (cost != INT_MAX)      
        {
            caminho_encontrado = 1; 
        }       
    }
    return make_pair(caminho_encontrado, cost);
}

int mapeamento(Hash h, string sequence, int k)
{
    int posicao = 0, length, aux, caminho_encontrado = 1;
    string resposta = "";
    vector<int> posicoesValidas;
    list<string> kmer_lista_aux;
    Hash dbg_gap(k);
    pair<int, int> status;
    int cost_aux = sequence.length();


    for (int i = 0; i < sequence.length() - (k - 1); i++)
    {
        if (h.contains(sequence.substr(i,k)))
        {
            posicoesValidas.push_back(i);
        }
    }

    cout << "Qtd. Anchros " << posicoesValidas.size() << endl;

    if (posicoesValidas.size() > 0)
    {
        if (utils.typeGraph == 0)
        {
            h.dbgToTraditionalSequenceGraph(0, 0);
            h.dbgToTraditionalSequenceGraph(1, 0);
        } else {
            h.dbgToSimplifiedSequenceGraph(0);
            h.dbgToSimplifiedSequenceGraph(1);
        }

        status = headMapping(posicoesValidas[0], h, sequence, k);
        if (status.second >= INT_MAX)
        {   
            int extra_cost = sequence.length() - posicoesValidas[0];
            return status.second + extra_cost;
        }

        // cout << "Head path " << status.first << endl;

        cost_aux = status.second;

        status = internalMapping(posicoesValidas, h, sequence, k);
        if (status.second >= INT_MAX)
        {
            int extra_cost = sequence.length() - status.first;
            return cost_aux + status.second + extra_cost;
        }

        // cout << "Internal path " << status.first << endl;

        cost_aux += status.second;

        status = tailMapping(posicoesValidas[posicoesValidas.size() - 1], h, sequence, k);
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
        
    Hash h(utils.k);   
    ifstream file(utils.nameSequenceArchive);
    h.populateGraph(utils.nameArchive, false);       

    while(getline(file, line))
    {
        //utils.readSequence(utils.nameSequenceArchive);  
        getline(file, line);
        cout << "Size L.Read " << line.size() << endl;
        utils.sequence = line;
        // mapeamento
        auto retorno = mapeamento(h, utils.sequence, utils.k); 
        cout << retorno << endl << endl;
    }     
    return 0;
}
