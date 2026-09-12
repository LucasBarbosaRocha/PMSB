# variante2 — Mapeamento de Sequências em Grafo de De Bruijn (com mudanças no grafo)

Esta pasta implementa a **variante 2** do PMSB (Capítulo 6 da tese,
`../Tese_doutorado.pdf`): em vez de editar a sequência `s` pra aproximá-la do
grafo `G` (variante 1, `../variante1/`), aqui é o **grafo** que é editado —
trocando o rótulo (k-mer) de alguns vértices — até que exista um percurso em
`G` que induz exatamente `s`. Custo de uma edição = soma das distâncias de
Hamming entre cada rótulo trocado e seu k-mer novo.

## O que muda em relação a Gibney et al.

Gibney et al. (2022) definiram uma versão desse problema em que a
**estrutura do grafo (arcos) é fixa** — só se pode trocar rótulos de vértices
sem que isso afete quais vértices são vizinhos de quais. Nessa versão, eles
provaram que o problema é **NP-completo**.

Este trabalho relaxa exatamente essa restrição: como as arestas de um grafo
de De Bruijn são **implícitas** (dois vértices são vizinhos se o sufixo do
k-mer de um é o prefixo do k-mer do outro), trocar o rótulo de um vértice
naturalmente **induz novos arcos e remove outros** — a topologia muda junto
com o rótulo, de graça, sem precisar de nenhuma contabilidade extra. Com essa
liberdade a mais, o problema deixa de ser NP-completo: a tese mostra que ele
se reduz a um **emparelhamento bipartido de custo mínimo** (k-mers da
sequência de um lado, vértices do grafo do outro, custo = distância de
Hamming), resolvido em tempo polinomial pelo **algoritmo Húngaro**. Esse
algoritmo é chamado na tese de **BMTC** (De Bruijn sequence Mapping Tool with
graph Changes).

Em resumo: **mesma pergunta, restrição diferente** (arcos fixos vs. arcos
induzidos) — e essa única diferença muda a classe de complexidade de
NP-completo para polinomial.

## Programas

| Binário/arquivo | O que é | Implementa o BMTC (Húngaro)? |
|---|---|---|
| `btl/btl.cpp` | "graph **L**abel change" | **Sim** — é a implementação real do Cap. 6: monta a matriz de custo (Hamming) entre os k-mers distintos de `s` e os vértices de `G`, roda o algoritmo Húngaro (`btl/Hungarian.cpp`, de terceiros) e aplica a edição de custo mínimo ao grafo. |
| `btn.cpp` / `btn_b.cpp` (com Bifrost) | "graph **N**ode change" | Não — é um protótipo anterior que só insere k-mers ausentes como vértices novos, sem custo nem emparelhamento. Mantido por referência histórica; para reproduzir o algoritmo da tese use `btl/`. |

## Compilando e rodando `btl` (recomendado)

```
cd variante2/btl
g++ -O2 -std=c++17 -c btl.cpp -o main.o
g++ -O2 -std=c++17 -c Hungarian.cpp -o hung.o
g++ -o btl main.o hung.o
```
(equivalente ao `make` já presente em `btl/makefile`, mas com `-std=c++17`)

```
./btl -s sequence2.fasta -g grafoT.fasta -k 3 -t 0
```

* `-s` arquivo FASTA com a sequência a ser mapeada
* `-g` arquivo FASTA com as sequências que formam o grafo de De Bruijn
* `-k` inteiro > 1, comprimento do k-mer
* `-t` tipo do grafo (mantido por compatibilidade com `myUtils`, não afeta o BMTC)

### Saída

* `Kmers no grafo: <n>` — tamanho do grafo (`|Gk|`)
* `Comprimento da sequência: <n>` — tamanho de `s`
* `Quantidade de Kmers na sequência presente no grafo (antes do mapeamento): <m>`
  — `|k(s)|`, quantidade de k-mers **distintos** de `s`
* Se `|Gk| < |k(s)|` (não é possível cobrir todos os k-mers de `s`): mensagem
  de erro, nenhum emparelhamento é tentado
* `Custo: <n>` — custo do emparelhamento ótimo (= custo da edição ótima do
  grafo, pelo Teorema 4 da tese)
* `Quantidade de Kmers na sequência presente no grafo (após o mapeamento): <n>`
  — quantas **janelas** (posições, não k-mers distintos) de `s` passam a
  existir no grafo depois da edição. É uma contagem diferente da anterior
  (uma conta k-mers distintos, essa conta todas as janelas/posições), então
  não é incomum os dois números divergirem mesmo num mapeamento correto —
  o valor esperado após um mapeamento bem-sucedido é o total de janelas de
  `s` (`comprimento - k + 1`), mostrando que agora *todas* estão cobertas.

## Compilando e rodando `btn`/`btn_b` (protótipo, não o algoritmo da tese)

```
cd variante2
g++ -O2 -std=c++17 -o btn btn.cpp
./btn -s <sequencia.fasta> -g <grafo.fasta> -k <k> -t 0
```

`btn_b.cpp` depende de Bifrost (`utils/myBifrost/`), que não está presente
nesta pasta — precisa copiar/linkar a mesma estrutura usada em
`../variante1/bmt/bmt_h1_b.cpp` para compilá-lo.

## Dados de exemplo

`btl/grafoT.fasta` + `btl/sequence2.fasta` (ou `sequence3.fasta`) são o par
de teste mais realista disponível na pasta; `btl/grafo.fasta`/`grafo2.fasta`
são grafos minúsculos (3 sequências) que servem só pra testar o caminho de
erro ("`Gk` pequeno demais pra `s`").
