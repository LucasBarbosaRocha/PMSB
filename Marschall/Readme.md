# Algoritmo de Rautiainen and Marschall

## Dado um grafo de sequências G e uma sequência s, o algoritmo devolve um percurso p em G tal que a sequência induzida por p é a mais parecida com s


## Grafo de De Bruijn

## Grafo de sequências simples

## Trasformando um grafo de De Bruijn em um grafo de sequências simples

### Ideia 1 (função dbgToSequenceGraph_1)

### Ideia 2 (função dbgToSequenceGraph_2)

## Mapeando uma sequência s no grafo de sequência simples

## Começando

### Pré-requisitos

* GCC (>= 4.0)
* MPI

## Rautiainen and Marschall

O arquivo do marschall.cpp possui três parâmetros de entrada: um grafo de De Bruijn (um arquivo listando sequências), uma sequência e um inteiro k >= 0.

### Entrada

Um arquivo txt com as sequências do Grafo de De Bruijn, uma sequência s e um inteiro positivo k.

**Exemplo de grafo de De Bruijn:**

```
Sequence 1
Sequence 2
Sequence 3
```
A sequência s representa a sequência mapeada no grafo de De Bruijn e o inteiro _k_ representa o comprimento do _k_-mer.

### Saída
A sequência induzida pelo caminho no grafo de De Bruijn tal que a diferença entre s e a sequência induzida é a menor possível.

### Execução
```
./marschall -s sequence -g graph.txt -k 3

## Resultados
|   |   |   |   |   |
|---|---|---|---|---|
|   |   |   |   |   |
|   |   |   |   |   |
|   |   |   |   |   |
