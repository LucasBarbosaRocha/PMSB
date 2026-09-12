# webui — painel local pra compilar e executar o PMSB

Site local (roda só na sua máquina, sem instalar nada — usa só a biblioteca
padrão do Python 3) pra: escolher uma ferramenta, compilar com um clique,
subir os arquivos de grafo e sequência, apertar "Executar" e ver o
resultado na tela. Toda execução é salva automaticamente em
`webui/results/<timestamp>_<ferramenta>.txt`. Grafo sempre
tradicional (`-t 0`) e mapeamento sempre por passeio (`-p 0`) — fixo, não é
escolha na interface.

## Rodando

```
python3 webui/server.py
```

Abre em **http://localhost:8765**.

## Ferramentas disponíveis

Ordem exibida no dropdown do site:

| # | Label no site | Sigla na tese | Ferramenta real | Pasta | Custo reportado? |
|---|---|---|---|---|---|
| 1 | Variante 1 (algoritmo exato) — BSMT: custo e mapeamento | BSMT | `bmt` | `variante1/bmt/` | sim (`Custo: N`) |
| 2 | Variante 1 (algoritmo exato) — BSMTd: apenas custo | BSMTd | `bmt_only_cost/bmt` | `variante1/bmt_only_cost/` | sim (`Custo: N`) — precisa da Bifrost compilada |
| 3 | Variante 1 (Heurística 1) — BSMTh1: custo e mapeamento | BSMTh1 | `bmt_h1` | `variante1/bmt/` | sim (`Custo: N`) |
| 4 | Variante 1 (Heurística 1 com Bifrost) — BSMTh1b: custo e mapeamento | BSMTh1b | `bmt_h1_b` | `variante1/bmt/` | sim (`Custo: N`) — precisa da Bifrost compilada |
| 5 | Variante 1 (Heurística 2) — BSMTh2: custo e mapeamento | BSMTh2 | `bmt_h2` | `variante1/bmt/` | sim (`Custo: N`) |
| 6 | Variante 1 (Heurística 3) — BSMTh3: custo e mapeamento | BSMTh3 | `bmt_h3` | `variante1/bmt/` | sim (`Custo: N`) |
| 7 | Variante 1 (Heurística 1) — BSMTh1d: apenas custo | BSMTh1d | `bmt_only_cost/bmt_h1` | `variante1/bmt_only_cost/` | sim (`Custo: N`) — precisa da Bifrost compilada |
| 8 | Variante 1 (Heurística 1 com Bifrost) — BSMTh1db: apenas custo | BSMTh1db | `bmt_only_cost/bmt_h1_b` | `variante1/bmt_only_cost/` | sim (`Custo: N`) — precisa da Bifrost compilada |
| 9 | Variante 2 — BMTC (emparelhamento Húngaro) | BMTC | `btl` | `variante2/btl/` | sim (`Custo: N`) |
| 10 | Auxiliar — Pesb | — (utilitário, não é um dos métodos formais) | `pesb` | `variante1/bmt/` | não tem custo, mostra k-mers ausentes (`Vamos inserir: N`) |

Saída de cada ferramenta padronizada: `Comprimento` (tamanho da sequência),
`Quantidade de Âncoras` (quando a ferramenta usa sementes), a sequência
mapeada rotulada como "Sequência mapeada (com as alterações aplicadas)"
(só nas ferramentas de custo+mapeamento) e `Custo: N` por último. `btl`
(BMTC) tem rótulos próprios (`Kmers no grafo`, `Comprimento da sequência`,
`Quantidade de Kmers na sequência presente no grafo` antes/depois do
mapeamento, `Custo: N`) porque mede coisas diferentes das ferramentas da
variante 1.

`variante2/btn` (protótipo "graph Node change", sem custo/emparelhamento —
não implementa o BMTC da tese) foi removido do catálogo do site por decisão
do usuário; o arquivo `btn.cpp` continua no repositório, só não aparece mais
como opção aqui.

As ferramentas "precisa da Bifrost compilada" compilam e rodam pelo site,
mas só depois da lib estar buildada localmente (não vem pronta no
repositório):

```
cd variante1/utils/myBifrost/bifrost
mkdir -p build && cd build
cmake .. && make -j$(nproc)
```

Se `libbifrost.a` não existir, o botão Compilar avisa isso e não tenta
compilar — sem quebrar o resto do site. Todas essas ferramentas compilam com
`-march=native`, que é **obrigatório, não cosmético**: a lib Bifrost usa
AVX2 internamente (headers como `BlockedBloomFilter.hpp` têm campos de
struct condicionados a `#if defined(__AVX2__)`), e compilar sem essa mesma
flag gera um binário que linka normalmente mas **crasha em runtime** por
corrupção de memória silenciosa (ODR violation, não erro de compilação).
Detalhes em `variante1/bmt/README.md`.

**De fora, de propósito**: `btn_b` (variante2) — a pasta `variante2/utils/`
nem tem o wrapper `myBifrost/` copiado, então não cabe no fluxo de um clique
sem antes replicar essa estrutura.

## O que cada botão faz

* **Compilar** — roda o(s) comando(s) `g++` daquela ferramenta na pasta
  correta e mostra o log de compilação (ou o erro, se falhar).
* **Executar** — sobe os dois arquivos pra `<pasta-da-ferramenta>/uploads/`,
  roda o binário já compilado com `-s -g -k -t 0 -p 0` (quando a ferramenta
  usa `-t`/`-p`), mostra a saída completa na tela, destaca o custo (quando a
  ferramenta reporta um — pega a última ocorrência de `Custo: N` na saída) e
  salva tudo — comando executado, parâmetros e saída — em `webui/results/`.

## Limitações conhecidas

* Servidor single-user, pensado pra uso local (não exponha `0.0.0.0` numa
  rede sem revisar segurança — ele roda binários compilados a partir dos
  arquivos `.cpp` do repositório, mas os *dados* enviados por upload não são
  sanitizados além do nome do arquivo).
* Usa o módulo `cgi` da biblioteca padrão para parsear o upload
  (`multipart/form-data`); esse módulo está descontinuado e será removido no
  Python 3.13 — se isso quebrar no futuro, trocar `cgi.parse_multipart` por
  um parser manual de multipart ou por uma dependência como `python-multipart`.
* `bmt_only_cost/bmt_h1_b` tem um leak pequeno e conhecido (~476 bytes por
  execução, ponteiro interno da Bifrost nunca liberado) — não afeta o
  resultado, só acumula memória em execuções muito longas/repetidas.
