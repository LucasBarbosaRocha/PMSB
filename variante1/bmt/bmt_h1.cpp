#include <fstream>
#include <iostream>
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"
#include <queue>
#include <bits/stdc++.h>

constexpr double kGapLengthFactor = 0.5;
constexpr int kUseHeuristic = 1;

MyUtils utils;
Marschall m;

list<int> nodes;
pair<list<string>, string> mapping;
pair<vector<pair<int,string>>, int> retorno;

int returnLengthGap(int k, int gap)
{
    if (gap > k)
        return gap * kGapLengthFactor;
    return k;
}

tuple<int, string, int> headMapping(Hash &dbg_gap, int firstPosition, Hash &h, const string &sequence, int k)
{
    int length, final, final_aux, aux, caminho_encontrado = 0, custo = 0;
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

        // buscando os kmers do gap
        final_aux = h.verifyLabelExistisAndReturnVerticeIndex(kmer_cauda);

        if (final_aux == -1)
        {
            caminho_encontrado = -1; resposta = "";
            return make_tuple(caminho_encontrado, resposta, custo);
        }

        final_aux += (k-1);

        nodes = h.sequenceGraphReverse.bfs(final_aux, length);
        dbg_gap.insertKmer(kmer_cauda, 0);        
        h.insertKmersByNodes(nodes, dbg_gap);

        // TODO CONFERIR os nós do do grafo simplificado
        if (utils.typeGraph == 0)
            dbg_gap.dbgToTraditionalSequenceGraph(0, !kUseHeuristic); 
        else
            dbg_gap.dbgToSimplifiedSequenceGraph(0); 

        string sequence_aux = sequence.substr(0, k + firstPosition);
        m.buildMultilayerGraph(dbg_gap.sequenceGraph, sequence_aux);        
        inicial_final = dbg_gap.vertice_inicial_final(kmer_cauda, kmer_cauda);
        final = (inicial_final.second + 2) + ((sequence_aux.size() - 1) * dbg_gap.sequenceGraph.getV()) + (sequence_aux.size() - 1) + (k-1);
        
        retorno = m.dijkstra(m.m_sequenceGraph, m.getInitialNode(), m.getEndNode());

        if (retorno.second == INT_MAX)
        {
            caminho_encontrado = -1;
            return make_tuple(caminho_encontrado, resposta, custo);
        } else {

            if (utils.typeGraph == 0)
                mapping = m.showTraditionalMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);
            else
                mapping = m.showSimplifiedMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);

            resposta = mapping.second.substr(0,firstPosition) + sequence.substr(firstPosition, k);
            custo = retorno.second;
        }

        // liberando memória
        retorno.first.clear();
        mapping.first.clear();
        mapping.second.clear();
    }
    return make_tuple(caminho_encontrado, resposta, custo);
}

tuple<int, string, int> internalMapping(Hash &dbg_gap, int pos, const vector<int> &positions, Hash &h, const string &sequence, int k)
{
    int caminho_encontrado = positions.size(), length, aux, posicao = positions[0], inicial_aux, final_aux, dif, custo = 0;
    string resposta = "", kmer_cabeca = "", kmer_cauda = "", sequence_aux = "";
    pair<int, int> inicial_final;
    // cout << "INTERNAL " << endl;

    for(int i = pos; i < positions.size() - 1; i++)
    {     
        caminho_encontrado = i + 1;
        dif = positions[i + 1] - positions[i];

        if (dif > k)
        {
            length = returnLengthGap(k, dif), aux = 0;
            kmer_cabeca = sequence.substr(positions[i], k);
            kmer_cauda = sequence.substr(positions[i+1], k);
            // buscando os kmers do gap
            inicial_aux = h.verifyLabelExistisAndReturnVerticeIndex(kmer_cabeca);
            final_aux = h.verifyLabelExistisAndReturnVerticeIndex(kmer_cauda);

            if (final_aux == -1 || inicial_aux == -1)
                return make_tuple(caminho_encontrado, resposta, custo);
            
            final_aux += (k-1);

            dbg_gap.deleteHash();

            dbg_gap.insertKmer(kmer_cabeca, 0);
            dbg_gap.insertKmer(kmer_cauda, 0);

            nodes = h.sequenceGraph.bfs(inicial_aux, length);
            h.insertKmersByNodes(nodes, dbg_gap);

            nodes = h.sequenceGraphReverse.bfs(final_aux, length);
            h.insertKmersByNodes(nodes, dbg_gap);

            if (utils.typeGraph == 0)
                dbg_gap.dbgToTraditionalSequenceGraph(0, !kUseHeuristic);        
            else
                dbg_gap.dbgToSimplifiedSequenceGraph(0);

            sequence_aux = sequence.substr(positions[i], k + dif);
            m.buildMultilayerGraph(dbg_gap.sequenceGraph, sequence_aux);
            inicial_final = dbg_gap.vertice_inicial_final(kmer_cabeca, kmer_cauda);
            inicial_aux = inicial_final.first + 2;
            final_aux = (inicial_final.second + 2 + (k-1)) + ((sequence_aux.size() - 1) * dbg_gap.sequenceGraph.getV()) + (sequence_aux.size() - 1);
            retorno = m.dijkstra(m.m_sequenceGraph, inicial_aux, final_aux);

            if (retorno.second == INT_MAX)
            {
                int posicao_atual = i;
                i = positions.size() + 1;
                if (positions[posicao_atual] == 0)
                    resposta = resposta + kmer_cabeca;
                else
                    resposta = resposta + kmer_cabeca[k-1];
                return make_tuple(caminho_encontrado, resposta, custo);
            } else
            {
                if (utils.typeGraph == 0)
                    mapping = m.showTraditionalMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);
                else
                    mapping = m.showSimplifiedMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);

                if (mapping.second.size() > k)
                {
                    // resposta = resposta + sequence.substr(positions[i], k);
                    resposta = resposta + mapping.second.substr(k, mapping.second.length() - k);
                    posicao = 1;
                    custo += retorno.second;
                }
            }

            // limpando memória
            dbg_gap.sequenceGraph.deleteGraph();
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
    return make_tuple(caminho_encontrado, resposta, custo);
}

tuple<int, string, int> tailMapping(Hash &dbg_gap, int lastPosition, Hash &h, const string &sequence, int k)
{
    int length, aux = 0, caminho_encontrado = -1, inicial_aux, dif, inicial, custo = 0;
    string resposta = "", kmer_cabeca = "", sequence_aux = "";
    pair<int, int> inicial_final;

    if (sequence.length() - lastPosition >= k)
    {
        // cout << "TAIL " << endl;
        length = returnLengthGap(k, sequence.length() - lastPosition);
        kmer_cabeca = sequence.substr(lastPosition, k);
        dbg_gap.insertKmer(kmer_cabeca, 0);   

        // buscando os kmers do gap
        inicial_aux = h.verifyLabelExistisAndReturnVerticeIndex(kmer_cabeca);
    
        if (inicial_aux == -1)
        {
            caminho_encontrado = -1; resposta = "";
            for (int trash_it = lastPosition + k; trash_it < sequence.length(); trash_it++)
                resposta = resposta + '-';
            return make_tuple(caminho_encontrado, resposta, custo);
        }
        dbg_gap.deleteHash();
        dbg_gap.insertKmer(kmer_cabeca, 0);
        nodes = h.sequenceGraph.bfs(inicial_aux, length);
        h.insertKmersByNodes(nodes, dbg_gap);

        if (utils.typeGraph == 0)
            dbg_gap.dbgToTraditionalSequenceGraph(0, !kUseHeuristic);
        else
            dbg_gap.dbgToSimplifiedSequenceGraph(0);

        dif = sequence.length() - lastPosition;    
        sequence_aux = sequence.substr(lastPosition, k + dif);
        m.buildMultilayerGraph(dbg_gap.sequenceGraph, sequence_aux);
        inicial_final = dbg_gap.vertice_inicial_final(kmer_cabeca, kmer_cabeca);
        inicial = inicial_final.first + 2;
        retorno = m.dijkstra(m.m_sequenceGraph, m.getInitialNode(), m.getEndNode());

        if (utils.typeGraph == 0)
            mapping = m.showTraditionalMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);
        else
            mapping = m.showSimplifiedMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);
    
        if (mapping.second.length() > k)
        {
            resposta = resposta + mapping.second.substr(k, mapping.second.length() - k);
            // cout << "tail: " << mapping.second.substr(k, mapping.second.length() - k);
            caminho_encontrado = 1;
            custo = retorno.second;
        } else {
            //for (int trash_it = lastPosition + k; trash_it < sequence.length(); trash_it++)
            //    resposta = resposta + '-';
        }
        retorno.first.clear();
        mapping.first.clear();
        mapping.second.clear();
    }
    return make_tuple(caminho_encontrado, resposta, custo);
}

pair<string, int> mapeamento(Hash &h, const string &sequence, int k)
{
    int posicao = 0, length, aux, caminho_encontrado = 1, custo_total = 0;
    string resposta = "";
    vector<int> posicoesValidas;
    Hash dbg_gap(k);
    tuple<int, string, int> status;


    for (int i = 0; i < sequence.length() - (k - 1); i++)
    {
        if (h.contains(sequence.substr(i,k)))
        {
            posicoesValidas.push_back(i);
        }
    }

    cout << "Quantidade de Âncoras: " << posicoesValidas.size() << endl;

    if (posicoesValidas.size() > 0)
    {
        if (utils.typeGraph == 0)
        {
            h.dbgToTraditionalSequenceGraph(0, kUseHeuristic);
            h.dbgToTraditionalSequenceGraph(1, kUseHeuristic);
        } else {
            h.dbgToSimplifiedSequenceGraph(0);
            h.dbgToSimplifiedSequenceGraph(1);
        }

        //h.sequenceGraph.printGraph();
        //h.sequenceGraphReverse.printGraph();

        string resp_temp = "";
        int custo_temp = 0;
        status = headMapping(dbg_gap, posicoesValidas[0], h, sequence, k);

        if (get<0>(status) != -1){
            resposta = get<1>(status);
            resp_temp = get<1>(status);
            custo_temp = get<2>(status);
        }

        int limite = posicoesValidas.size();

        int pos_internal = 0;
        string aux = resposta;
        int aux_custo = custo_temp;
        while (get<0>(status) < limite)
        {
            status = internalMapping(dbg_gap, pos_internal, posicoesValidas, h, sequence, k);
            if (get<0>(status) < limite)
            {
                aux = aux + get<1>(status);
                aux_custo = aux_custo + get<2>(status);
                if (resp_temp.size() < aux.size())
                {
                    resp_temp = aux;
                    custo_temp = aux_custo;
                }
            }
            pos_internal = get<0>(status);
            aux = "";
            aux_custo = 0;
        }

        resposta = resp_temp;
        custo_total = custo_temp;
        status = tailMapping(dbg_gap, posicoesValidas[posicoesValidas.size() - 1], h, sequence, k);
        if (get<0>(status) != -1) {
            resposta = resposta + get<1>(status);
            custo_total += get<2>(status);
        }
    }

    if (resposta.size() > 0)
        return make_pair(resposta, custo_total);
    else
        return make_pair("Sequência nao pode ser mapeada", custo_total);
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
        cout << "Comprimento: " << line.size() << endl;
        utils.sequence = line;
        // mapeamento
        auto retorno = mapeamento(h, utils.sequence, utils.k);
        cout << "Sequência mapeada (com as alterações aplicadas): " << retorno.first << endl;
        cout << "Custo: " << retorno.second << endl << endl;
    }
    return 0;
}
