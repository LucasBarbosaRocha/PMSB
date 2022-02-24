#include <iostream>
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"
#include <queue>
#include <bits/stdc++.h>

#define x 0.5
MyUtils utils;
Marschall m;

pair<int, string> headMapping(int firstPosition, Hash h, string sequence, int k)
{
    Hash dbg_gap(k);
    int length, aux, caminho_encontrado = 1;
    list<string> kmer_lista_aux;
    list<int> nodes;
    string resposta;
    Marschall m;
    pair<list<string>, string>  mapeamento;

    if (firstPosition != 0)
    {
        cout << "HEAD " << endl;
        length = k * x, aux = 0;
        string kmer_cauda = sequence.substr(firstPosition, k);

        // buscando os kmers do gap
        int final_aux = h.verifyLabelExistisAndReturnVerticeIndex(kmer_cauda);
    
        if (final_aux == -1)
        {
            caminho_encontrado = 0; resposta = "";
            return make_pair(caminho_encontrado, resposta);
        }

        nodes = h.sequenceGraphReverse.bfs(final_aux, length);
        dbg_gap.insertKmer(kmer_cauda);        
        h.insertKmersByNodes(nodes, dbg_gap);

        // TODO CONFERIR os nós do do grafo simplificado
        if (utils.typeGraph == 0)
            dbg_gap.dbgToTraditionalSequenceGraph(0); 
        else
            dbg_gap.dbgToSimplifiedSequenceGraph(0); 

        string sequence_aux = sequence.substr(0, k + firstPosition);
        m.buildMultilayerGraph(dbg_gap.sequenceGraph, sequence_aux);        
        auto inicial_final = dbg_gap.vertice_inicial_final(kmer_cauda, kmer_cauda);
        int final = (inicial_final.second + 2) + ((sequence_aux.size() - 1) * dbg_gap.sequenceGraph.getV()) + (sequence_aux.size() - 1);
        
        auto retorno = m.dijkstra(m.m_sequenceGraph, m.getInitialNode(), m.getEndNode());

        if (retorno.second == INT_MAX)  
        {
            caminho_encontrado = 0;
            return make_pair(caminho_encontrado, resposta);
        }

         if (utils.typeGraph == 0)
            mapeamento = m.showTraditionalMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);
        else
            mapeamento = m.showSimplifiedMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);

        resposta = mapeamento.second.substr(0,firstPosition) + sequence.substr(firstPosition, k);    

        // liberando memória
        retorno.first.clear();
        mapeamento.second.clear();
    }
    return make_pair(caminho_encontrado, resposta);
}

pair<int, string> internalMapping(vector<int> positions, Hash h, string sequence, int k)
{
    int caminho_encontrado = 0, length, aux, posicao = positions[0];
    string resposta;
    list<int> nodes;
    Hash dbg_gap(k);
    Marschall m;
    pair<list<string>, string> mapeamento;
    cout << "INTERNAL " << endl;
    for(int i = 0; i < positions.size() - 1; i++)
    {     
        caminho_encontrado = 1;

        int dif = positions[i + 1] - positions[i]; 
        if (dif > 1)
        {
            length = k * x, aux = 0;
            string kmer_cabeca = sequence.substr(positions[i], k);
            string kmer_cauda = sequence.substr(positions[i+1], k);
            // buscando os kmers do gap
            int inicial_aux = h.verifyLabelExistisAndReturnVerticeIndex(kmer_cabeca);
            int final_aux = h.verifyLabelExistisAndReturnVerticeIndex(kmer_cauda);
        
            if (final_aux == -1 || inicial_aux == -1)
            {
                caminho_encontrado = 0; resposta = "";
                return make_pair(caminho_encontrado, resposta);
            }
            dbg_gap.deleteHash();
            dbg_gap.insertKmer(kmer_cabeca);
            dbg_gap.insertKmer(kmer_cauda);
            nodes = h.sequenceGraph.bfs(inicial_aux, length);
            h.insertKmersByNodes(nodes, dbg_gap);
            nodes = h.sequenceGraphReverse.bfs(final_aux, length);
            h.insertKmersByNodes(nodes, dbg_gap);

            cout << "dbg_gap " << dbg_gap.getQtdKmers() << endl;

            if (utils.typeGraph == 0)
                dbg_gap.dbgToTraditionalSequenceGraph(0);        
            else
                dbg_gap.dbgToSimplifiedSequenceGraph(0);

            string sequence_aux = sequence.substr(positions[i], k + dif);
            m.buildMultilayerGraph(dbg_gap.sequenceGraph, sequence_aux);
            auto inicial_final = dbg_gap.vertice_inicial_final(kmer_cabeca, kmer_cauda);
            int inicial = inicial_final.first + 2;
            int final = (inicial_final.second + 1) + ((sequence_aux.size() - 1) * dbg_gap.sequenceGraph.getV()) + (sequence_aux.size() - 1);
            auto retorno = m.dijkstra(m.m_sequenceGraph, m.getInitialNode(), m.getEndNode());

            if (utils.typeGraph == 0)
                mapeamento = m.showTraditionalMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);
            else
                mapeamento = m.showSimplifiedMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);

            if (retorno.second == INT_MAX)
            {
                i = positions.size() + 1;
                caminho_encontrado = 0;
                if (positions[i] == 0)
                    resposta = resposta + kmer_cabeca;
                else
                    resposta = resposta + kmer_cabeca[k-1];
                break;
            } else
            if (mapeamento.second.size() > k)      
            {
                resposta = resposta + mapeamento.second.substr(k, mapeamento.second.length() - (2 * k));   
                resposta = resposta + sequence.substr(positions[i+1], k);                 
                posicao = 1;
            }

            // limpando memória
            dbg_gap.sequenceGraph.deleteGraph();
            mapeamento.second.clear();
            retorno.first.clear();
        }else{
            if (posicao == 0)
                resposta = resposta + sequence.substr(positions[i], k);
            else
                resposta = resposta + sequence.substr(positions[i] + (k-1), 1);
            posicao++;
        }
    }
    return make_pair(caminho_encontrado, resposta);
}

pair<int, string> tailMapping(int lastPosition, Hash h, string sequence, int k)
{
    int length = k * x, aux = 0, caminho_encontrado = 0;
    Hash dbg_gap(k);
    Marschall m;
    list<int> nodes;
    string resposta = "";
    pair<list<string>, string> mapeamento;

    if (sequence.length() - lastPosition >= k)
    {
        cout << "TAIL " << endl;

        string kmer_cabeca = sequence.substr(lastPosition, k);
        dbg_gap.insertKmer(kmer_cabeca);   

        // buscando os kmers do gap
        int inicial_aux = h.verifyLabelExistisAndReturnVerticeIndex(kmer_cabeca);
    
        if (inicial_aux == -1)
        {
            caminho_encontrado = 0; resposta = "";
            return make_pair(caminho_encontrado, resposta);
        }
        dbg_gap.deleteHash();
        dbg_gap.insertKmer(kmer_cabeca);
        nodes = h.sequenceGraph.bfs(inicial_aux, length);
        h.insertKmersByNodes(nodes, dbg_gap);

        if (utils.typeGraph == 0)
            dbg_gap.dbgToTraditionalSequenceGraph(0);
        else
            dbg_gap.dbgToSimplifiedSequenceGraph(0);

        int dif = sequence.length() - lastPosition;    
        string sequence_aux = sequence.substr(lastPosition, k + dif);
        m.buildMultilayerGraph(dbg_gap.sequenceGraph, sequence_aux);
        auto inicial_final = dbg_gap.vertice_inicial_final(kmer_cabeca, kmer_cabeca);
        int inicial = inicial_final.first + 2;
        auto retorno = m.dijkstra(m.m_sequenceGraph, m.getInitialNode(), m.getEndNode());

        if (utils.typeGraph == 0)
            mapeamento = m.showTraditionalMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);
        else
            mapeamento = m.showSimplifiedMapping(retorno.first, dbg_gap, dbg_gap.sequenceGraph);
    
        if (mapeamento.second.length() >= k)      
        {
            resposta = resposta + mapeamento.second.substr(k, mapeamento.second.length() - k); 
            caminho_encontrado = 1; 
        }       
    }
    return make_pair(caminho_encontrado, resposta);
}

string mapeamento(Hash h, string sequence, int k)
{
    int posicao = 0, length, aux, caminho_encontrado = 1;
    string resposta = "";
    vector<int> posicoesValidas;
    list<string> kmer_lista_aux;
    Marschall m;
    Hash dbg_gap(k);
    pair<int, string> status;


    for (int i = 0; i < sequence.length() - (k - 1); i++)
    {
        if (h.contains(sequence.substr(i,k)))
        {
            posicoesValidas.push_back(i);
        }
    }

    if (posicoesValidas.size() > 0)
    {
        if (utils.typeGraph == 0)
        {
            h.dbgToTraditionalSequenceGraph(0);
            h.dbgToTraditionalSequenceGraph(1);
        } else {
            h.dbgToSimplifiedSequenceGraph(0);
            h.dbgToSimplifiedSequenceGraph(1);
        }

        status = headMapping(posicoesValidas[0], h, sequence, k);
        if (status.first == 0)
            return status.second;

        cout << "Head path " << status.first << endl;

        resposta = status.second;

        status = internalMapping(posicoesValidas, h, sequence, k);
        if (status.first == 0)
            return resposta + status.second;

        cout << "Internal path " << status.first << endl;

        resposta += status.second;

        status = tailMapping(posicoesValidas[posicoesValidas.size() - 1], h, sequence, k);
        if (status.first == 0)
            return resposta + status.second;

        
        cout << "Tail path " << status.first << endl;
    } 

    if (resposta.size() > 0)
        return resposta;
    else
        return "Sequência nao pode ser mapeada";
}

int main(int argc, char *argv[])
{
    if(utils.verifyData(argc, argv) == 1)
        exit (0); 
        
    Hash h(utils.k);   
    utils.readSequence(utils.nameSequenceArchive);  
    h.populateGraph(utils.nameArchive, false);       
    // mapeamento
    auto retorno = mapeamento(h, utils.sequence, utils.k); 
    cout << retorno << endl;
    return 0;
}
