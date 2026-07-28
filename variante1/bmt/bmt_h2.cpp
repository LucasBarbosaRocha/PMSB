#include <fstream>
#include <iostream>
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"
#include <queue>
#include <bits/stdc++.h>

MyUtils utils;
Marschall m;

// Distância de edição (Levenshtein) entre duas sequências, com custo 1 para
// substituição, inserção e deleção — a mesma métrica usada pelo bmt/BSMT.
// Usada para reportar um custo real e comparável entre a heurística e o
// método exato (bmt_h2/bmt_h3 não usam Dijkstra, então não têm um custo
// acumulado pronto como bmt_h1).
int distanciaEdicao(const string &a, const string &b)
{
    int n = a.size(), m = b.size();
    vector<vector<int>> dist(n + 1, vector<int>(m + 1));
    for (int i = 0; i <= n; i++) dist[i][0] = i;
    for (int j = 0; j <= m; j++) dist[0][j] = j;
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            dist[i][j] = min({
                dist[i-1][j] + 1,
                dist[i][j-1] + 1,
                dist[i-1][j-1] + (a[i-1] != b[j-1] ? 1 : 0)
            });
    return dist[n][m];
}

// ESTENDE (Pseudocódigo 2 da tese): dada uma sequência q, um grafo de De
// Bruijn (h) e k, devolve a sequência induzida pelo maior prefixo de q cujos
// k-mers consecutivos existem todos em h. Se o prefixo inteiro de q é válido,
// o resultado tem o mesmo tamanho de q.
string estende(Hash &h, const string &q, int k)
{
    string mapping = "";
    if ((int)q.size() < k)
        return mapping;

    int limite = (int)q.size() - k; // último início de janela válido em q
    for (int i = 0; i <= limite; i++)
    {
        string kmer = q.substr(i, k);
        if (!h.contains(kmer))
            break;
        if (mapping.empty())
            mapping = kmer;
        else
            mapping += kmer.substr(k - 1, 1);
    }
    return mapping;
}

// Heurística 2 (BSMT_h2) da tese: ancora sementes (k-mers de s presentes no
// grafo), e para cada gap entre duas sementes consecutivas tenta, para cada
// base do alfabeto, substituir o último caractere da semente esquerda e
// estender via ESTENDE até alcançar a semente direita. Se nenhuma base
// alcança a semente direita, guarda o melhor resultado parcial encontrado e
// recomeça a partir do próximo gap (mantendo o melhor mapeamento global já
// obtido, seguindo o critério de maior comprimento).
pair<string, int> seed_and_extend(Hash &h, const string &sequence, int k)
{
    const string bases = "ACGT";
    vector<int> positions;

    for (int i = 0; i + k <= (int)sequence.size(); i++)
        if (h.contains(sequence.substr(i, k)))
            positions.push_back(i);

    cout << "Qtd. Anchros " << positions.size() << endl;

    if (positions.empty())
        return make_pair("Sequência nao pode ser mapeada", (int)sequence.length());

    string p_best = "";
    string p_temp = sequence.substr(positions[0], k);

    for (size_t idx = 0; idx + 1 < positions.size(); idx++)
    {
        int a = positions[idx];
        int aLinha = positions[idx + 1];
        int dif = aLinha - a;

        if (dif <= k)
        {
            // Sementes sobrepostas/adjacentes: nenhuma busca é necessária,
            // os caracteres entre elas já são garantidos pela sobreposição
            // dos dois k-mers válidos. Faltam exatamente "dif" caracteres
            // (de a+k até aLinha+k-1) para p_temp alcançar a nova semente.
            p_temp += sequence.substr(a + k, dif);
            continue;
        }

        string q = sequence.substr(a, (aLinha + k) - a);
        bool alcancei = false;
        string melhorLocal = "";

        for (char c : bases)
        {
            string qMod = q;
            qMod[k - 1] = c;
            string p = estende(h, qMod, k);

            if (p.size() == qMod.size())
            {
                alcancei = true;
                melhorLocal = p;
                break;
            }
            if (p.size() >= melhorLocal.size())
                melhorLocal = p;
        }

        if (alcancei)
        {
            // melhorLocal cobre [a, aLinha+k); os primeiros k caracteres já
            // estão representados no fim de p_temp, então só o restante é novo.
            p_temp += melhorLocal.substr(k);
        }
        else
        {
            if (p_temp.size() >= p_best.size())
                p_best = p_temp;
            p_temp = melhorLocal;
        }
    }

    if (p_temp.size() >= p_best.size())
        p_best = p_temp;

    if (p_best.size() > 0)
        return make_pair(p_best, distanciaEdicao(sequence, p_best));
    else
        return make_pair("Sequência nao pode ser mapeada", (int)sequence.length());
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
        getline(file, line);
        cout << "Size L.Read " << line.size() << endl;
        transform(line.begin(), line.end(), line.begin(), ::toupper);
        utils.sequence = line;
        auto retorno = seed_and_extend(h, utils.sequence, utils.k);
        cout << retorno.first << endl;
        cout << retorno.second << endl;
    }
    return 0;
}
