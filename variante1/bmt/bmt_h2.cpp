#include <fstream>
#include <iostream>
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"
#include "../utils/heuristicsCommon.cpp"
#include <queue>
#include <bits/stdc++.h>

MyUtils utils;
Marschall m;

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

    cout << "Quantidade de Âncoras: " << positions.size() << endl;

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
        cout << "Comprimento: " << line.size() << endl;
        transform(line.begin(), line.end(), line.begin(), ::toupper);
        utils.sequence = line;
        auto retorno = seed_and_extend(h, utils.sequence, utils.k);
        cout << "Sequência mapeada (com as alterações aplicadas): " << retorno.first << endl;
        cout << "Custo: " << retorno.second << endl << endl;
    }
    return 0;
}
