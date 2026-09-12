/*
    Funções compartilhadas pelas heurísticas BSMT_h2 e BSMT_h3.
    Requer que "marschall.cpp" (ou qualquer header que declare a classe Hash)
    já tenha sido incluído antes deste arquivo.
*/
#pragma once
#include <string>
#include <vector>

// Distância de edição (Levenshtein) entre duas sequências, com custo 1 para
// substituição, inserção e deleção — a mesma métrica usada pelo bmt/BSMT.
// Usada para reportar um custo real e comparável entre a heurística e o
// método exato (bmt_h2/bmt_h3 não usam Dijkstra, então não têm um custo
// acumulado pronto como bmt_h1).
int distanciaEdicao(const std::string &a, const std::string &b)
{
    int n = a.size(), m = b.size();
    std::vector<std::vector<int>> dist(n + 1, std::vector<int>(m + 1));
    for (int i = 0; i <= n; i++) dist[i][0] = i;
    for (int j = 0; j <= m; j++) dist[0][j] = j;
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            dist[i][j] = std::min({
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
std::string estende(Hash &h, const std::string &q, int k)
{
    std::string mapping = "";
    if ((int)q.size() < k)
        return mapping;

    int limite = (int)q.size() - k; // último início de janela válido em q
    for (int i = 0; i <= limite; i++)
    {
        std::string kmer = q.substr(i, k);
        if (!h.contains(kmer))
            break;
        if (mapping.empty())
            mapping = kmer;
        else
            mapping += kmer.substr(k - 1, 1);
    }
    return mapping;
}
