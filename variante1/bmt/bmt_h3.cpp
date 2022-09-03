#include <fstream>
#include <iostream>
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"
#include <queue>
#include <bits/stdc++.h>

MyUtils utils;
Marschall m;

#define x 0.3
#define heuristc 1
unordered_map<int, string> mapp; 
unordered_map<int, int> vald; 

pair<int, string> findKmerInSequenceByValidPositions(string sequence, int pos, int k)
{
    int i = 0, pos_aux = pos;
    string kmer = "";
    cout << "validations " << ": ";
    while(pos_aux < sequence.size() && i < k)
    {
        cout << pos_aux << ": " << vald[pos_aux] << " ";
        if (vald[pos_aux] == 1)
        {
            kmer = kmer + sequence.substr(pos_aux, 1);
            i++;
        }
        pos_aux++;
    }
    cout << endl;
    return make_pair(pos_aux, kmer);
}

pair<int, string> findValidBaseInSequence(string sequence, int pos)
{
    int find = 1; string ch = "";
    while (find == 1 && pos < sequence.size())
    {
        if (vald[pos] == 1)
        {
            find = 0;
            ch = sequence.substr(pos, 1);
        } else
            pos++;
    }
    return make_pair(pos, ch);
}

void cleanMappAndVald(int pos, int k)
{       
    cout << "clean " << endl;                 
    int temp_pos = pos;
    for(int temp_k = 0; temp_k < k; temp_k++)
    {
        mapp[temp_pos] = '-';
        cout << vald[temp_pos] << " -> ";
        vald[temp_pos] = 0;
        cout << vald[temp_pos] << endl;
        temp_pos++;
    }
}

tuple<string, string, int> extend(Hash h, string kmer, string sequence, int pos, int k)
{
    string mapping = "", ch;
    pair<int, string> pos_base;
    if (h.contains(kmer))
    {
        mapping = kmer; 
        cout << "Comecando com pos " << pos << endl;
        pos_base = findValidBaseInSequence(sequence, pos);
        pos = pos_base.first; ch = pos_base.second;
        kmer = kmer.substr(1,k-1) + ch;

        pos++;
        while (1)
        {
            if (h.contains(kmer))
            {
                mapping += ch;

                pos_base = findValidBaseInSequence(sequence, pos);
                pos = pos_base.first; ch = pos_base.second;
                kmer = kmer.substr(1,k-1) + ch;

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
    cout << "pos " << pos << ": " << kmer_busca << endl;
    tuple <string, string, int> saida = extend(h, kmer_busca, sequence, pos, k);
    string retorno = get<0>(saida);
    return saida;
}

pair<string, int> seed_and_extend(Hash h, string sequence, int k)
{
    int posicao = 0, length, aux, caminho_encontrado = 1, errors = 0, faltaMapeamento = 0, pos = 0, pos_tmp = 0;
    string resposta = "", resposta_local = "", resposta_temp = "";
    vector<int> positions;
    list<string> kmer_lista_aux;
    Marschall m;
    Hash dbg_gap(k);
    pair<int, string> status;
    string bases = {'A', 'C', 'G', 'T'}, kmer = "";

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

    for (int i = 0; i < limit; i++)
    {
        mapp[i] = sequence.substr(i, 1);
        vald[i] = 1;
    }

    while (pos < limit)
    {
        // kmer = sequence.substr(pos, k);
        cout << "Pos " << pos << endl;
        status = findKmerInSequenceByValidPositions(sequence, pos, k);
        pos = status.first; kmer = status.second;

        cout << "kmer " << kmer << endl;

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
                mapp[pos+k] = 'A';
            }

            saida = extend_from_kmer_with_new_base(h, kmer_aux, "C", sequence, pos, k);
            retorno = get<0>(saida);

            if (retorno.size() > tam_atual)
            {
                resposta_atual = retorno;
                tam_atual = retorno.size();
                pos_aux = get<2>(saida);
                kmer_aux = get<1>(saida);
                mapp[pos+k] = 'C';
            }

            saida = extend_from_kmer_with_new_base(h, kmer_aux, "G", sequence, pos, k);
            retorno = get<0>(saida);

            if (retorno.size() > tam_atual)
            {
                resposta_atual = retorno;
                tam_atual = retorno.size();
                pos_aux = get<2>(saida);
                kmer_aux = get<1>(saida);
                mapp[pos+k] = 'G';
            }

            saida = extend_from_kmer_with_new_base(h, kmer_aux, "T", sequence, pos, k);
            retorno = get<0>(saida);
        
            if (retorno.size() > tam_atual)
            {
                resposta_atual = retorno;
                tam_atual = retorno.size();
                pos_aux = get<2>(saida);
                kmer_aux = get<1>(saida);
                mapp[pos+k] = 'T';
            }

            if (pos == pos_aux) // nao consegui estender, vou apagar k caracteres para trás
            {
                cout << "hora de apagar " << pos << " " << k << endl;
                if (sequence.size() > k)
                {
                    if (pos - k > 0)                    
                        cleanMappAndVald(pos - k, k);                   
                    else
                        cleanMappAndVald(pos, k);

                    if (pos - k > 0)
                        pos -= k; 
                    limit -= k;
                    //if (resposta.size() > k && resposta.size() - (k - 1) > 0)
                    //    resposta.erase(resposta.size() - k - 1, k - 1);
                } else
                    break;
            } else {
                resposta = resposta + resposta_atual;  
                pos = pos_aux;
            }
        } else {
            if (pos == 0)
                resposta = kmer;
            else {
                resposta += kmer.substr(k-1,1);
            }
            pos++; 
        }    
    }
    for (int i = 0; i < mapp.size(); i++)
    {
        cout << vald[i] << ": " << mapp[i] << endl;
    }

    resposta = "";
    for (int i = 0; i < mapp.size(); i++)
    {
        resposta = resposta + mapp[i];
    }

    if (resposta.size() > 0)
        return make_pair(resposta, errors);
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
