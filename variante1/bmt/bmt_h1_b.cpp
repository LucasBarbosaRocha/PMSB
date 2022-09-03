#include <fstream>
#include <iostream>
#include "../utils/myUtils.cpp"
#include "../utils/myBifrost/my_bifrost.cpp"
#include <queue>
#include <bits/stdc++.h>

#define x 0.5
#define heuristc 1

MyUtils utils;

list<int> nodes;
pair<list<string>, string> mapping;
pair<vector<pair<int,string>>, int> retorno;

int returnLengthGap(int k, int gap)
{
    if (gap > k)
        return gap * x;
    return k;
}

pair<int, string> headMapping(my_Bifrost dbg_gap, int firstPosition, my_Bifrost h, string sequence, int k)
{
    int length, final, final_aux, aux, caminho_encontrado = 0;
    list<string> kmer_lista_aux;
    pair<int, int> inicial_final;
    string resposta = "";

    if (firstPosition != 0)
    {
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
        auto grafo_multicada = dbg_gap.buildMultilayerGraph(dbg_gap.sequenceGraph, sequence_aux);
        auto retorno = dbg_gap.dijkstra(grafo_multicada.first, 0, grafo_multicada.second);
        mapping = dbg_gap.buildMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);   

        if (retorno.second == INT_MAX)  
        {
            caminho_encontrado = -1;
            return make_pair(caminho_encontrado, resposta);
        } else {
            resposta = mapping.second.substr(0,firstPosition) + sequence.substr(firstPosition, k);    
        }
        // liberando memória
        dbg_gap.clear();
        retorno.first.clear();
        mapping.first.clear();
        mapping.second.clear();
    }
    return make_pair(caminho_encontrado, resposta);
}

pair<int, string> internalMapping(my_Bifrost dbg_gap, int pos, vector<int> positions, my_Bifrost h, string sequence, int k)
{
    int caminho_encontrado = positions.size(), length, aux, posicao = positions[0], inicial_aux, final_aux, dif;
    string resposta = "", kmer_cabeca = "", kmer_cauda = "", sequence_aux = "";
    pair<int, int> inicial_final;
    // cout << "INTERNAL " << endl;

    for(int i = pos; i < positions.size() - 1; i++)
    {     
        caminho_encontrado = i + 1;
        dif = positions[i + 1] - positions[i];

        if (dif > k)
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
            auto grafo_multicada = dbg_gap.buildMultilayerGraph(dbg_gap.sequenceGraph, sequence_aux);

            auto inicial_final = dbg_gap.vertice_inicial_final(kmer_cabeca, kmer_cauda);
            auto inicial_teste = inicial_final.first + 2;
            auto final_teste = (inicial_final.second + 2 + (k-1)) + ((sequence_aux.size() - 1) * dbg_gap.sequenceGraph.getV()) + (sequence_aux.size() - 1);
            inicial_aux = 0;
            final_aux = (grafo_multicada.second + 2 + (k-1)) + ((sequence_aux.size() - 1) * dbg_gap.sequenceGraph.getV()) + (sequence_aux.size() - 1);


            //auto retorno = dbg_gap.dijkstra(grafo_multicada.first, inicial_teste, grafo_multicada.second);           
            auto retorno = dbg_gap.dijkstra(grafo_multicada.first, inicial_teste, final_teste);           
            mapping = dbg_gap.buildMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);   
            // cout << "retorno " << retorno.second << endl;
            // cout << "mapping " << mapping.second << endl;
            if (retorno.second == INT_MAX)
            {
                i = positions.size() + 1;
                if (positions[i] == 0)
                    resposta = resposta + kmer_cabeca;
                else
                    resposta = resposta + kmer_cabeca[k-1]; 
                return make_pair(caminho_encontrado, resposta);
            } else
            {
                if (mapping.second.size() > k)      
                {
                    // resposta = resposta + sequence.substr(positions[i], k);
                    resposta = resposta + mapping.second.substr(k, mapping.second.length() - k);   
                    posicao = 1;
                }
            }

            // limpando memória
            dbg_gap.sequenceGraph.deleteGraph();
            dbg_gap.clear();
            mapping.first.clear();
            mapping.second.clear();
            retorno.first.clear();
        }else{
            //cout << positions[i] << " " << positions[i+1] << endl; 
            if (posicao == 0)
                resposta = resposta + sequence.substr(positions[i], k);
            else {
                resposta = resposta + sequence.substr(positions[i] + (k-1), 1);
                if (positions[i + 1] == sequence.length() - k)
                    resposta = resposta + sequence.substr(positions[i + 1] + (k-1), 1);
            }
            //cout << resposta << endl;
            posicao++;
        }
    }
    return make_pair(caminho_encontrado, resposta);
}

pair<int, string> tailMapping(my_Bifrost dbg_gap, int lastPosition, my_Bifrost h, string sequence, int k)
{
    int length, aux = 0, caminho_encontrado = -1, inicial_aux, dif, inicial;
    string resposta = "", kmer_cabeca = "", sequence_aux = "";
    pair<int, int> inicial_final;

    if (sequence.length() - lastPosition >= k)
    {
        // cout << "TAIL " << endl;
        length = returnLengthGap(k, sequence.length() - lastPosition);
        kmer_cabeca = sequence.substr(lastPosition, k);

        dbg_gap.insertSequence(kmer_cabeca);
        h.traverseByDistance(dbg_gap, kmer_cabeca, length, 1);

        dbg_gap.convertDbgToSequenceGraph();
        dif = sequence.length() - lastPosition;    
        sequence_aux = sequence.substr(lastPosition, k + dif);

        auto grafo_multicada = dbg_gap.buildMultilayerGraph(dbg_gap.sequenceGraph, sequence_aux);

        auto retorno = dbg_gap.dijkstra(grafo_multicada.first, 0, grafo_multicada.second);
        mapping = dbg_gap.buildMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);   
    
        if (mapping.second.length() > k)      
        {
            resposta = resposta + mapping.second.substr(k, mapping.second.length() - k); 
            // cout << "tail: " << mapping.second.substr(k, mapping.second.length() - k);
            caminho_encontrado = 1; 
        } else {
            //for (int trash_it = lastPosition + k; trash_it < sequence.length(); trash_it++)
            //    resposta = resposta + '-';    
        }            
        
        dbg_gap.clear();
        retorno.first.clear();
        mapping.first.clear();
        mapping.second.clear();
    }
    return make_pair(caminho_encontrado, resposta);
}

string mapeamento(my_Bifrost h, string sequence, int k)
{
    int posicao = 0, length, aux, caminho_encontrado = 1;
    string resposta = "";
    vector<int> posicoesValidas;
    my_Bifrost dbg_gap(k);
    pair<int, string> status;


    posicoesValidas = h.findAnchors(sequence);
    cout << "Qtd. Anchros " << posicoesValidas.size() << endl;
    /* for (auto it : posicoesValidas)
        cout << it << " ";
    cout << endl; */
    if (posicoesValidas.size() > 0)
    {
        if (utils.typeGraph == 0)
        {

            // h.dbgToTraditionalSequenceGraph(0, heuristc);
            // h.dbgToTraditionalSequenceGraph(1, heuristc);
        } else {
            //h.dbgToSimplifiedSequenceGraph(0);
            //h.dbgToSimplifiedSequenceGraph(1);
        }

        //h.sequenceGraph.printGraph();
        //h.sequenceGraphReverse.printGraph();

        string resp_temp = "";
        status = headMapping(dbg_gap, posicoesValidas[0], h, sequence, k);
        if (status.first != -1){
            resposta = status.second;
            resp_temp = status.second;
        }

        int limite = posicoesValidas.size();

        int pos_head = 0, pos_internal = 0, pos_internal_escolhida = 0;
        string aux = resposta;
        while (status.first < limite)
        {
            status = internalMapping(dbg_gap, pos_internal, posicoesValidas, h, sequence, k);
            if (status.first < limite) 
            {
                aux = aux + status.second;
                if (resp_temp.size() < aux.size())
                {
                    resp_temp = aux;
                    pos_internal_escolhida = pos_internal;
                }
            }  
            pos_internal = status.first;
            aux = "";
        }

        resposta = resp_temp;
        status = tailMapping(dbg_gap, posicoesValidas[posicoesValidas.size() - 1], h, sequence, k);
        if (status.first != -1) {    
            resposta = resposta + status.second;
        } 
    } 

    if (resposta.size() > 0)
        return resposta;
    else
        return "Sequência nao pode ser mapeada";
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
        cout << "Kmers: " << bf.size() << endl;
        cout << "Size L.Read " << line.size() << endl;
        utils.sequence = line;
        // mapeamento
        auto retorno = mapeamento(bf, utils.sequence, utils.k); 
        cout << retorno << endl << endl;
    }     
    return 0;
}
