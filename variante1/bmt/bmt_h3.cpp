#include <fstream>
#include <iostream>
#include "../utils/marschall.cpp"
#include "../utils/myUtils.cpp"
#include "../utils/heuristicsCommon.cpp"
#include <queue>
#include <bits/stdc++.h>

MyUtils utils;
Marschall m;

// Passo da Heurística 2 aplicado a um único gap: tenta as 4 bases no último
// caractere do k-mer inicial de q e estende. Devolve o melhor resultado e se
// a extensão alcançou o fim de q (ou seja, a semente da direita).
pair<string, bool> bsmtH2Gap(Hash &h, const string &q, int k)
{
    const string bases = "ACGT";
    string melhor = "";
    bool alcancei = false;

    for (char c : bases)
    {
        string qMod = q;
        qMod[k - 1] = c;
        string p = estende(h, qMod, k);

        if (p.size() == qMod.size())
        {
            alcancei = true;
            melhor = p;
            break;
        }
        if (p.size() >= melhor.size())
            melhor = p;
    }
    return make_pair(melhor, alcancei);
}

// Heurística 3 (BSMT_h3) da tese: para cada gap entre sementes consecutivas,
// tenta primeiro o mesmo passo da Heurística 2. Se não alcançar a semente da
// direita, tenta remover k-1 caracteres logo após o primeiro caractere da
// semente esquerda (uma exclusão) e estender novamente a partir daí,
// mantendo o resultado que cobrir mais.
pair<string, int> seed_and_extend(Hash &h, const string &sequence, int k)
{
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
            p_temp += sequence.substr(a + k, dif);
            continue;
        }

        string q = sequence.substr(a, (aLinha + k) - a);
        auto [p_local, alcancei] = bsmtH2Gap(h, q, k);

        if (!alcancei && (int)q.size() > k)
        {
            // Não alcançamos a semente da direita: tenta remover k-1
            // caracteres logo após o 1o caractere da semente esquerda
            // (uma exclusão) e estender de novo a partir daí.
            string qCompacta = q.substr(0, 1) + q.substr(k);
            string p2 = estende(h, qCompacta, k);
            if (p2.size() >= p_local.size())
            {
                p_local = p2;
                alcancei = (p2.size() == qCompacta.size());
            }
        }

        if (alcancei)
        {
            p_temp += p_local.substr(k);
        }
        else
        {
            if (p_temp.size() >= p_best.size())
                p_best = p_temp;
            else if (p_local.size() >= p_best.size())
                p_best = p_local;
            p_temp = "";
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
