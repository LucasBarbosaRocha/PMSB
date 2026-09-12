# bmt — Mapeamento de Sequências em Grafo de De Bruijn

Esta pasta contém as implementações do PMSB (Problema do Mapeamento de
Sequências em grafo de De Bruijn) da variante 1 da tese ("O Problema do
Mapeamento de Sequências em Grafos de De Bruijn", Lucas Barbosa Rocha,
UFMS/FACOM, 2024 — `../../Tese_doutorado.pdf`). Dados uma sequência `s` e um
grafo de De Bruijn `G`, o objetivo é encontrar o percurso em `G` cuja
sequência induzida seja a mais parecida possível com `s`. Se você só precisa
do **custo** do mapeamento (não da sequência mapeada em si), veja a variante
mais leve em `../bmt_only_cost/`.

Os programas desta pasta implementam o algoritmo exato (Capítulo 5.1 da
tese, **BSMT**) e as três heurísticas propostas (Capítulo 5.2, **BSMT_h1**,
**BSMT_h2**, **BSMT_h3**). Veja o `README.md` da pasta `variante1/` para o
contexto geral do trabalho e a publicação associada.

## Programas

| Binário    | Arquivo         | Método (nome na tese)                                                          |
|------------|------------------|----------------------------------------------------------------------------------|
| `bmt`      | `bmt.cpp`        | **BSMT** — algoritmo exato via Dijkstra no grafo de multicamadas (referência)     |
| `bmt_h1`   | `bmt_h1.cpp`     | **Heurística 1 (BSMT_h1)** — ancora sementes e resolve cada gap com Dijkstra num subgrafo local |
| `bmt_h1_b` | `bmt_h1_b.cpp`   | Igual ao `bmt_h1`, usando a biblioteca Bifrost para construir o grafo de De Bruijn |
| `bmt_h2`   | `bmt_h2.cpp`     | **Heurística 2 (BSMT_h2)** — seed-and-extend com substituição                    |
| `bmt_h3`   | `bmt_h3.cpp`     | **Heurística 3 (BSMT_h3)** — seed-and-extend com substituição e exclusão         |
| `pesb`     | `pesb.cpp`       | Utilitário: conta quantos k-mers da sequência não existem no grafo               |

`bmt` é o método exato (mais lento, usado como referência de qualidade —
**nenhuma heurística pode encontrar um custo menor que o de `bmt`** para a
mesma entrada, já que `bmt` é ótimo). `bmt_h1`, `bmt_h1_b`, `bmt_h2` e
`bmt_h3` são heurísticas mais rápidas, com qualidade de mapeamento variável
dependendo da entrada — inclusive podendo devolver mapeamentos **parciais**
(mais curtos que a sequência de entrada, ou com gaps `-`) quando não
conseguem cobrir a sequência inteira. Isso é um comportamento esperado e
documentado na tese, não um bug.

## Compilando

```
./build.sh
```

Compila `bmt`, `bmt_h1`, `bmt_h2`, `bmt_h3` e `pesb` (dependem apenas dos
arquivos em `../utils/`) e coloca os binários em `bmt/bin/`.

`bmt_h1_b` depende da biblioteca [Bifrost](https://github.com/pmelsted/bifrost)
(código-fonte em `../utils/myBifrost/bifrost/`), que **não vem pré-compilada**
neste repositório. Compile-a primeiro:

```
cd ../utils/myBifrost/bifrost
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

Depois rode `./build.sh` de novo — ele detecta `libbifrost.a` e linka
`bmt_h1_b` automaticamente. Se a lib ainda não existir, ele avisa e segue sem
interromper o build dos demais programas.

**Atenção com `-march=native`:** o `CMakeLists.txt` da Bifrost compila a lib
com `COMPILATION_ARCH=native` (ativa AVX2 nesta máquina). Vários headers da
Bifrost (ex. `BlockedBloomFilter.hpp`) têm campos de struct condicionados a
`#if defined(__AVX2__)` — se `bmt_h1_b.cpp` for compilado **sem** essa mesma
flag, o tamanho dessas structs diverge entre o `.cpp` e a lib já compilada
(violação de ODR), e o binário **linka normalmente mas crasha em runtime**
com corrupção de memória silenciosa (não dá erro de compilação nem de link —
só aparece rodando, e o sintoma parece aleatório: SEGV numa estrutura sem
relação nenhuma com o bug real). Por isso `build.sh` já usa `-march=native`
ao compilar `bmt_h1_b` — se compilar manualmente, use a mesma flag.

## Entrada

Todo programa recebe os mesmos parâmetros de linha de comando:

```
<binário> -s nameSequenceArchive -g nameGraphArchive -k kmerSize -t typeSequenceGraph -p pathOrWalk
```

* `-s` arquivo com a sequência a ser mapeada
* `-g` arquivo com as sequências que formam o grafo de De Bruijn
* `-k` inteiro > 1, comprimento do k-mer
* `-t` tipo do grafo de sequências: `0` tradicional, `1` simplificado
* `-p` `0` para passeio (walk) ou `1` para caminho (path) durante o mapeamento

Ambos os arquivos (`-s` e `-g`) usam formato FASTA simples — um cabeçalho
começando com `>` seguido de uma linha com a sequência de bases (A, C, G, T).
O arquivo do grafo (`-g`) pode ter várias sequências; cada uma contribui
k-mers para o grafo de De Bruijn. O arquivo da sequência (`-s`) pode ter
várias sequências também; cada uma é mapeada, uma de cada vez.

Exemplo (`examples/grafo.fasta`):
```
>ref1
ACGTACGTTGCA
>ref2
TTGCATGCACGT
```

Exemplo (`examples/sequence.fasta`):
```
>read1
ACGTCCGTTGCA
```

## Saída

A saída varia um pouco entre os programas, mas em geral cada um imprime:

* `Comprimento: <n>`: tamanho da sequência lida
* `Quantidade de Âncoras: <n>` (heurísticas): quantidade de k-mers da
  sequência que batem exatamente com algum k-mer do grafo — as **sementes**
  usadas para guiar a busca (chamado de *seeds* / conjunto `A` na tese)
* `Sequência mapeada (com as alterações aplicadas): <seq>` — o resultado do
  algoritmo, ou seja, a sequência de entrada já com as edições
  (substituição/inserção/deleção) aplicadas para se encaixar no grafo. Pode
  ser mais curta que a entrada (mapeamento parcial) ou conter `-` no lugar
  de uma base (gap) quando a heurística não consegue cobrir a sequência
  inteira.
* `Custo: <n>`: custo do mapeamento, na mesma métrica de edição usada pelo
  algoritmo exato (substituição/inserção/deleção = custo 1 cada). Em `bmt`
  é sempre o custo **ótimo**. Em `bmt_h1` é a soma dos custos de Dijkstra de
  cada segmento resolvido. Em `bmt_h2`/`bmt_h3` (que não usam Dijkstra) é a
  distância de edição (Levenshtein) calculada entre a sequência de entrada
  e o resultado mapeado. Em todos os casos, o custo de uma heurística deve
  ser sempre **maior ou igual** ao custo de `bmt` para a mesma entrada —
  nunca menor, já que `bmt` é ótimo.

`pesb` tem saída diferente: mostra quantos k-mers da sequência **não**
existem no grafo (`Vamos inserir: <n>`), útil para diagnosticar rapidamente
quão "distante" uma sequência está do grafo antes de rodar um mapeamento
completo.

## Executando o exemplo

Com os binários já compilados (`./build.sh`), rode:

```
./run_examples.sh
```

Isso executa `bmt`, `bmt_h1`, `bmt_h2`, `bmt_h3` e `pesb` com os arquivos de
`examples/` (`k=3`) e mostra o comando e a saída de cada um. Se os binários
ainda não existirem, o script chama `./build.sh` automaticamente antes.

### Comandos básicos (equivalentes, chamados manualmente)

```
cd variante1/bmt
./build.sh

./bin/bmt    -s examples/sequence.fasta -g examples/grafo.fasta -k 3 -t 0 -p 0
./bin/bmt_h1 -s examples/sequence.fasta -g examples/grafo.fasta -k 3 -t 0 -p 0
./bin/bmt_h2 -s examples/sequence.fasta -g examples/grafo.fasta -k 3 -t 0 -p 0
./bin/bmt_h3 -s examples/sequence.fasta -g examples/grafo.fasta -k 3 -t 0 -p 0
./bin/pesb   -s examples/sequence.fasta -g examples/grafo.fasta -k 3 -t 0 -p 0
```

### Saída esperada no exemplo

A sequência de entrada (`ACGTCCGTTGCA`) tem uma base diferente de `ref1`
(`ACGTACGTTGCA`) na posição 4 (`C` em vez de `A`).

**`bmt`** (custo ótimo):
```
Comprimento: 12
Sequência mapeada (com as alterações aplicadas): ACGTACGTTGCA
Custo: 1
```
Encontra o caminho de custo mínimo no grafo e devolve exatamente `ref1`,
com custo 1 (uma substituição) — o resultado correto e ótimo.

**`bmt_h1`**:
```
Comprimento: 12
Quantidade de Âncoras: 7
Sequência mapeada (com as alterações aplicadas): ACG--TTGCA
Custo: 4
```
Mapeamento parcial (2 gaps) com custo 4 — pior que o ótimo (1), como
esperado de uma heurística, mas nunca melhor. A busca de `bmt_h1` é
limitada a um subgrafo local ao redor de cada gap (ver Seção 5.2.1 da
tese), então nem sempre encontra o caminho completo.

**`bmt_h2`**:
```
Comprimento: 12
Quantidade de Âncoras: 7
Sequência mapeada (com as alterações aplicadas): CGTTGCA
Custo: 5
```

**`bmt_h3`**:
```
Comprimento: 12
Quantidade de Âncoras: 7
Sequência mapeada (com as alterações aplicadas): TGCA
Custo: 8
```
`bmt_h2` e `bmt_h3` são heurísticas *seed-and-extend* mais simples (sem
Dijkstra): a partir de cada semente, tentam estender caractere a caractere
substituindo (e, no caso de `bmt_h3`, também excluindo) bases até perder o
alinhamento com o grafo, devolvendo o melhor trecho contíguo encontrado.
Por isso o resultado pode ser mais curto que a entrada — não é uma falha,
é o "melhor trecho" que essas heurísticas conseguiram ancorar no grafo
para este exemplo pequeno. O custo (5 e 8, respectivamente) é a distância
de edição entre a entrada e esse resultado parcial — maior que o custo
ótimo de `bmt` (1), como esperado de uma heurística.

## Dados maiores

Para testar com dados reais/maiores, veja `../Oficial/archives/` (exemplos
de `kmers.txt`/`sequence.txt`) ou gere seus próprios arquivos FASTA seguindo
o formato descrito acima. Tenha em mente que os métodos heurísticos
(`bmt_h1`/`bmt_h1_b`) reconstroem sub-grafos a cada gap encontrado, então o
tempo de execução cresce com a quantidade de gaps na sequência, não só com
o tamanho do grafo.
