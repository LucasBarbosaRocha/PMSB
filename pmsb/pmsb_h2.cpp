#include <fstream>
#include <iostream>
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"
#include <queue>
#include <bits/stdc++.h>

#define x 0.3
#define heuristc 1
MyUtils utils;
Marschall m;

tuple<string, string, int> extend(Hash h, string kmer, string sequence, int pos, int k)
{
    string mapping = "";
    if (h.contains(kmer))
    {
        mapping = kmer;
        kmer = kmer.substr(1,k-1) + sequence.substr(pos,1);
        pos++;
        while (1)
        {
            if (h.contains(kmer))
            {
                mapping += sequence.substr(pos,1);
                kmer = kmer.substr(1,k-1) + sequence.substr(pos,1);
                // cout << "- procura " << kmer << endl;
                pos = pos + 1;
            } else{
                break;
            }
        }
    } 
    return make_tuple(mapping, kmer, pos);
}

tuple<string, string, int> extend_from_kmer_with_new_base(Hash h, string kmer, string base, string sequence, int pos, int k)
{
    string kmer_busca = kmer + base;
    tuple <string, string, int> saida = extend(h, kmer_busca, sequence, pos, k);
    string retorno = get<0>(saida);

    // cout << "Comeca com " << kmer_busca << " : " << retorno << " last_kmer " <<  get<1>(saida) << " final " << kmer_cauda << endl;

    return saida;
}

pair<string, int> seed_and_extend(Hash h, string sequence, int k)
{
    int posicao = 0, length, aux, caminho_encontrado = 1, errors = 0, faltaMapeamento = 0, pos = 0, pos_oficial = 0, index = 0;
    string resposta = "", resposta_local = "", resposta_oficial = "";
    vector<int> positions;
    list<string> kmer_lista_aux;
    Marschall m;
    Hash dbg_gap(k);
    pair<int, string> status;
    string bases = {'A', 'C', 'G', 'T'};

    for (int i = 0; i < sequence.length() - (k - 1); i++)
    {
        if (h.contains(sequence.substr(i,k)))
        {
            positions.push_back(i);
        }
    }

    cout << "Qtd. Anchros " << positions.size() << endl;
    resposta = "";

    int limit = sequence.size();

    pos = 0;
    if (positions.size() > 0)
        pos = positions[index];
    faltaMapeamento = 1;
    while (pos < limit)
    {

        string kmer = sequence.substr(pos, k);
        if (!h.contains(kmer))
        {
            string kmer_aux = kmer.substr(0, k-1);
            int tam_atual = 0, pos_aux = pos;

            tuple <string, string, int> saida = extend_from_kmer_with_new_base(h, kmer_aux, "A", sequence, pos, k);
            string retorno = get<0>(saida), resposta_atual;
            if (retorno.size() > tam_atual)
            {
                resposta_atual = retorno;
                tam_atual = retorno.size();
                pos_aux = get<2>(saida);
                kmer_aux = get<1>(saida);
            }

            saida = extend_from_kmer_with_new_base(h, kmer_aux, "C", sequence, pos, k);
            retorno = get<0>(saida);

            if (retorno.size() > tam_atual)
            {
                resposta_atual = retorno;
                tam_atual = retorno.size();
                pos_aux = get<2>(saida);
                kmer_aux = get<1>(saida);
            }

            saida = extend_from_kmer_with_new_base(h, kmer_aux, "G", sequence, pos, k);
            retorno = get<0>(saida);

            if (retorno.size() > tam_atual)
            {
                resposta_atual = retorno;
                tam_atual = retorno.size();
                pos_aux = get<2>(saida);
                kmer_aux = get<1>(saida);
            }

            saida = extend_from_kmer_with_new_base(h, kmer_aux, "T", sequence, pos, k);
            retorno = get<0>(saida);
        
            if (retorno.size() > tam_atual)
            {
                resposta_atual = retorno;
                tam_atual = retorno.size();
                pos_aux = get<2>(saida);
                kmer_aux = get<1>(saida);
            }
            
            if (pos == pos_aux) // nao consegui estender, vou apagar k caracteres para trás
            {
                if (resposta.size() > resposta_oficial.size())
                {
                    resposta_oficial = resposta;
                    pos_oficial = pos;
                    faltaMapeamento = 1;
                    resposta = "";
                }
                pos++;
            } else {
                resposta = resposta + resposta_atual;  
                pos = pos_aux;
            }
        } else {
            if (faltaMapeamento == 1){
                resposta = kmer; faltaMapeamento = 0;}
            else {
                resposta += kmer.substr(k-1,1);}
            pos++; 
        }    
    }

    if (resposta.size() > resposta_oficial.size())
        resposta_oficial = resposta;

    if (resposta_oficial.size() > 0)
        return make_pair(resposta_oficial, errors);
    else
        return make_pair("Sequência nao pode ser mapeada", sequence.length());
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
        transform(line.begin(), line.end(), line.begin(), ::toupper);
        utils.sequence = line;
        // mapeamento
        auto retorno = seed_and_extend(h, utils.sequence, utils.k); 
        cout << retorno.first << endl;
        cout << retorno.second << endl;
    }     
    return 0;
}
