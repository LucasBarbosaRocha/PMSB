# bmt_only_cost — versão "só custo" (BSMTd / BSMTh1d / BSMTh1db)

Esta pasta é a versão "só distância" dos programas de `../bmt/`: em vez de
reconstruir e devolver a sequência mapeada, cada programa devolve **apenas o
custo** do mapeamento (a distância de edição mínima entre `s` e o percurso
encontrado em `G`). Na tese (`../../Tese_doutorado.pdf`) isso corresponde ao
**Algoritmo 7 (BSMTd)** e às heurísticas equivalentes **BSMTh1d**/**BSMTh1db**:
como não é preciso guardar o caminho para reconstruir a sequência, o algoritmo
usa uma estrutura de "duas camadas" (`../utils/marschallTwoLayers.cpp`, que
mantém só a camada anterior e a atual do grafo de multicamadas, não o grafo
inteiro) — ganhando memória em troca de não devolver o mapeamento em si. A
tese registra essa troca no Cap. 8: BSMTd consegue representar grafos maiores
que o BSMT completo, mas ainda assim pode levar várias horas em instâncias
grandes.

## Programas

| Binário      | Arquivo            | Equivale a (`../bmt/`) | Método (nome na tese) |
|--------------|---------------------|-------------------------|------------------------|
| `bmt`        | `bmt.cpp`           | `bmt`                   | **BSMTd** — só o custo do algoritmo exato |
| `bmt_h1`     | `bmt_h1.cpp`        | `bmt_h1`                | **BSMTh1d** — só o custo da Heurística 1 |
| `bmt_h1_b`   | `bmt_h1_b.cpp`      | `bmt_h1_b`              | **BSMTh1db** — igual, usando Bifrost |

Use esta pasta quando só precisa saber **quão perto** `s` está de `G`
(ex.: para comparar várias sequências rapidamente, ou para instâncias grandes
onde reconstruir o mapeamento completo seria caro). Use `../bmt/` quando
precisa da sequência mapeada em si.

## Compilando

Não existe `build.sh` nesta pasta ainda. `bmt.cpp` e `bmt_h1.cpp` não usam a
biblioteca [Bifrost](https://github.com/pmelsted/bifrost) diretamente, mas
`../utils/marschallTwoLayers.cpp` inclui `../utils/myBifrost/my_bifrost.cpp`
internamente — então o *link* dos três binários desta pasta (não só o
`bmt_h1_b`) exige a Bifrost já compilada.

Primeiro compile a Bifrost (uma vez só, se ainda não tiver feito):

```
cd variante1/utils/myBifrost/bifrost
mkdir -p build && cd build
cmake .. && make -j$(nproc)
```

Depois compile os três binários **com `-march=native`** (obrigatório, não
cosmético — veja a explicação em `../bmt/README.md`, seção "Atenção com
-march=native": a lib Bifrost usa AVX2 internamente, e compilar sem essa
flag gera um binário que linka mas crasha em runtime por corrupção de
memória silenciosa):

```
cd variante1/bmt_only_cost
mkdir -p bin

g++ -O2 -march=native -std=c++17 -o bin/bmt bmt.cpp \
  -I ../utils/myBifrost/bifrost/src \
  ../utils/myBifrost/bifrost/build/src/libbifrost.a -lpthread -lz

g++ -O2 -march=native -std=c++17 -o bin/bmt_h1 bmt_h1.cpp \
  -I ../utils/myBifrost/bifrost/src \
  ../utils/myBifrost/bifrost/build/src/libbifrost.a -lpthread -lz

g++ -O2 -march=native -std=c++17 -o bin/bmt_h1_b bmt_h1_b.cpp \
  -I ../utils/myBifrost/bifrost/src \
  ../utils/myBifrost/bifrost/build/src/libbifrost.a -lpthread -lz
```

## Entrada

Mesmo formato de `../bmt/` — grafo e sequência em FASTA simples:

```
<binário> -s nameSequenceArchive -g nameGraphArchive -k kmerSize -t typeSequenceGraph -p pathOrWalk
```

* `-s` arquivo com a sequência a ser mapeada
* `-g` arquivo com as sequências que formam o grafo de De Bruijn
* `-k` inteiro > 1, comprimento do k-mer
* `-t` tipo do grafo de sequências: `0` tradicional, `1` simplificado
* `-p` `0` para passeio (walk) ou `1` para caminho (path) — exigido pelo
  `myUtils.cpp` compartilhado, mesmo que estes binários não usem a distinção

Exemplos prontos nesta pasta: `grafo.fasta`, `sequence.fasta`.

```
./bin/bmt      -s sequence.fasta -g grafo.fasta -k 3 -t 0 -p 0
./bin/bmt_h1   -s sequence.fasta -g grafo.fasta -k 3 -t 0 -p 0
./bin/bmt_h1_b -s sequence.fasta -g grafo.fasta -k 3 -t 0 -p 0
```

## Saída

Bem mais enxuta que `../bmt/`:

* `Comprimento: <n>`: tamanho da sequência lida
* `Quantidade de Âncoras: <n>` (`bmt_h1`/`bmt_h1_b`): quantidade de sementes ancoradas
* `Custo: <n>` — o custo do mapeamento (distância de edição mínima entre `s`
  e o melhor percurso encontrado). Em `bmt` é sempre o custo **ótimo**; em
  `bmt_h1`/`bmt_h1_b` é uma estimativa heurística, sempre **maior ou igual**
  ao custo ótimo, refletindo que este é o produto do BSMTh1d/BSMTh1db na
  tese (Algoritmo 7: "devolve apenas a distância do mapeamento"). Nenhum dos
  três devolve a sequência mapeada em si.

Quando o mapeamento falha em cobrir parte da sequência (gap sem solução), o
custo devolvido é um limite superior de fallback (baseado no comprimento não
coberto), não `INT_MAX` — isso garante que o número impresso sempre é
utilizável para comparação, mesmo em mapeamentos parciais malsucedidos.
