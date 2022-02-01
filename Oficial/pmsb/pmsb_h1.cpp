#include <iostream>
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"
#include <queue>
#include <bits/stdc++.h>


string mapeamento(Hash h, string sequence, int k)
{
    int posicao = 0, length, aux, caminho_encontrado = 1;
    float x = 0.5;
    string resposta = "";
    vector<int> posicoesValidas;
    list<string> kmer_lista_aux;
    Marschall m;
    Hash dbg_gap(k);
    auto sequenceGraph = h.dbgToTraditionalSequenceGraph(0);       
    auto sequenceGraphReverse = h.dbgToTraditionalSequenceGraph(1);   

    for (int i = 0; i < sequence.length() - (k - 1); i++)
    {
        if (h.contains(sequence.substr(i,k)))
        {
            posicoesValidas.push_back(i);

        }
    }

    if (posicoesValidas.size() > 0)
    {
        if (posicoesValidas[0] != 0)
        {
            length = posicoesValidas[0] * x, aux = 0;
            string kmer_cauda = sequence.substr(posicoesValidas[0], k);

            dbg_gap.insertKmer(kmer_cauda);   
            // buscando os kmers do gap
            int final_aux = sequenceGraph.verifyLabelExistisAndReturnVerticeIndex(kmer_cauda);

            auto vertices = sequenceGraphReverse.bfs(final_aux, length);
            for (auto vertice : vertices.first)
            {
                dbg_gap.insertKmer(sequenceGraph.getBase(vertice.first));
                kmer_lista_aux.push_back(sequenceGraph.getBase(vertice.first));
            }

            auto sequence_graph_aux = dbg_gap.dbgToTraditionalSequenceGraph(0);    
        
            string sequence_aux = sequence.substr(0, k + posicoesValidas[0]);
            auto grafo_multicamada = m.buildMultilayerGraph(sequence_graph_aux, sequence_aux);
            auto inicial_final = h.vertice_inicial_final(kmer_cauda, kmer_cauda);
            int final = inicial_final.second + ((sequence_aux.size() - 1) * sequence_graph_aux.getV());
            auto retorno = m.dijkstra(grafo_multicamada, m.getInitialNode(), final);
            auto mapeamento = m.showTraditionalMapping(retorno.first, h, sequence_graph_aux);
            if (retorno.second == INT_MAX)
            {
                caminho_encontrado = 0;
            }else
                resposta = mapeamento.second.substr(0,posicoesValidas[0]) + sequence.substr(posicoesValidas[0], k); 
            posicao++;
        }


        if(caminho_encontrado == 1)
        {     
            for(int i = 0; i < posicoesValidas.size() - 1; i++)
            {     
                caminho_encontrado = 1;
                //cout << "(" << posicoesValidas[i] << "," << posicoesValidas[i + 1] << ") "; 
                int dif = posicoesValidas[i + 1] -  posicoesValidas[i]; 
                if (dif > 1)
                {

                    length = dif * x, aux = 0;
                    string kmer_cabeca = sequence.substr(posicoesValidas[i], k);
                    string kmer_cauda = sequence.substr(posicoesValidas[i+1], k);
                    dbg_gap.insertKmer(kmer_cabeca);
                    dbg_gap.insertKmer(kmer_cauda);
                    //buscando os kmers do gap
                    int inicial_aux = sequenceGraph.verifyLabelExistisAndReturnVerticeIndex(kmer_cabeca);
                    int final_aux = sequenceGraph.verifyLabelExistisAndReturnVerticeIndex(kmer_cauda);


                    auto vertices = sequenceGraph.bfs(inicial_aux, length);
                    cout << "aqui " << endl;

                    auto vertices2 = sequenceGraphReverse.bfs(final_aux, length);


                    for (auto vertice : vertices.first)
                    {
                        dbg_gap.insertKmer(sequenceGraph.getBase(vertice.first));
                    }
                    
                    for (auto vertice : vertices2.first)
                    {
                        dbg_gap.insertKmer(sequenceGraph.getBase(vertice.first));
                    }

                    auto sequence_graph_aux = dbg_gap.dbgToTraditionalSequenceGraph(0);        
                    string sequence_aux = sequence.substr(posicoesValidas[i], k + dif);
                    auto grafo_multicamada = m.buildMultilayerGraph(sequence_graph_aux, sequence_aux);
                    auto inicial_final = h.vertice_inicial_final(kmer_cabeca, kmer_cauda);
                    int inicial = inicial_final.first + 2;
                    int final = inicial_final.second + ((sequence_aux.size() - 1) * sequence_graph_aux.getV());
                    auto retorno = m.dijkstra(grafo_multicamada, inicial, final);
                    auto mapeamento = m.showTraditionalMapping(retorno.first, h, sequence_graph_aux);
                 

                    if (retorno.second == INT_MAX)
                    {
                        i = posicoesValidas.size() + 1;
                        caminho_encontrado = 0;
                        if (posicoesValidas[i] == 0)
                            resposta = resposta + kmer_cabeca;
                        else
                            resposta = resposta + kmer_cabeca[k-1];
                        break;
                    } else
                    if (mapeamento.second.size() >= k)      
                    {
                        resposta = resposta + mapeamento.second.substr(k, mapeamento.second.length() - (2 * k));   
                        //resposta = resposta + mapeamento.substr(mapeamento.length()- (k), k);   
                        resposta = resposta + sequence.substr(posicoesValidas[i+1], k);                 
                        posicao = 1;
                    }
                }else{
                    if (posicao == 0)
                        resposta = resposta + sequence.substr(posicoesValidas[i], k);
                    else
                        resposta = resposta + sequence.substr(posicoesValidas[i] + (k-1), 1);
                    posicao++;
                }
            }


            if (caminho_encontrado == 1 && posicoesValidas[posicoesValidas.size() - 1] != sequence.length() - k)
            {
                //cout << "(" << posicoesValidas[posicoesValidas.size() - 1] << ",-)" << endl;
                length = (sequence.length() - posicoesValidas[posicoesValidas.size() - 1]) * x, aux = 0;
                string kmer_cabeca = sequence.substr(posicoesValidas[posicoesValidas.size() - 1], k);
                dbg_gap.insertKmer(kmer_cabeca);   
                //kmers_mapeamento.push_back(kmer_cabeca);

                // buscando os kmers do gap
                int inicial_aux = sequenceGraph.verifyLabelExistisAndReturnVerticeIndex(kmer_cabeca);

                auto vertices = sequenceGraphReverse.bfs(inicial_aux, length);
                for (auto vertice : vertices.first)
                {
                    dbg_gap.insertKmer(sequenceGraph.getBase(vertice.first));
                    kmer_lista_aux.push_back(sequenceGraph.getBase(vertice.first));
                }
                auto sequence_graph_aux = dbg_gap.dbgToTraditionalSequenceGraph(0);
    
                int dif = sequence.length() - posicoesValidas[posicoesValidas.size() - 1];
                if (dif > 1)
                {    
                    string sequence_aux = sequence.substr(posicoesValidas[posicoesValidas.size() - 1], k + dif);
                    auto grafo_multicamada = m.buildMultilayerGraph(sequence_graph_aux, sequence_aux);
                    auto inicial_final = h.vertice_inicial_final(kmer_cabeca, kmer_cabeca);
                    int inicial = inicial_final.first + 2;
                    auto retorno = m.dijkstra(grafo_multicamada, inicial, m.getEndNode());
                    auto mapeamento = m.showTraditionalMapping(retorno.first, h, sequence_graph_aux);
                    if (mapeamento.second.length() >= k)      
                        resposta = resposta + mapeamento.second.substr(k, mapeamento.second.length() - k);  
                
                }
            }
        }
    } else
    {
        cout << "Nenhum kmer encontrado no grafo. A sequência vai ser toda mapeada." << endl;
        auto sequence_graph = h.dbgToTraditionalSequenceGraph(0);
        auto grafo_multicamada = m.buildMultilayerGraph(sequence_graph, sequence);
        auto retorno = m.dijkstra(grafo_multicamada, m.getInitialNode(), m.getEndNode());
        auto mapeamento = m.showTraditionalMapping(retorno.first, h, sequence_graph);
        return resposta;
    }

    if (resposta.size() > 0)
        return resposta;
    else
        return "caminho nao encontrado";
}

int main(int argc, char *argv[])
{
    MyUtils utils;
    if(utils.verifyData(argc, argv) == 1)
        exit (0);
    
    Hash h(utils.k);   
    utils.readSequence(utils.nameSequenceArchive);
    h.populateGraph(utils.nameArchive, false);
    // mapeamento
    cout << "qtd. kmers " << h.getQtdKmers() << endl;
    auto retorno = mapeamento(h, utils.sequence, utils.k); 
    cout << retorno << endl;
    return 0;
}
