# Trabalho de Doutorado 

# PMSB -- o Problema do Mapeamento de Sequências em grafo de De Bruijn (com mudanças na sequência).


## PMSB
Dadas uma sequência s e um grafo de De Bruijn G, encontrar o menor percurso p em G, tal que a sequência induzida por p seja a mais parecida com s.

### Entrada

Um arquivo txt com as sequências do Grafo de De Bruijn, uma sequência s e um inteiro positivo k.

**Exemplo de grafo de De Bruijn:**

```
Sequence 1
Sequence 2
Sequence 3
```
### Saída
A sequência induzida pelo caminho no grafo de De Bruijn tal que a diferença entre s e a sequência induzida é a menor possível.

### Execução
```
./bmt -s nameSequenceArchive -g nameGraphArchive -k kmerSize -t typeSequenceGraph (0 - traditional, 1 - simplified) -p (0 - walk, 1 - path)
```

* -s para a sequência a ser mapeada no grafo
* -g para o arquivo com o grafo de De Bruijn
* -k com o inteiro > 1 para o comprimento do k-mer
* -t para o tipo do grafo
* -p para usar passeio ou caminho durante o mapeamento


# Publicações
* Heuristics for the de Bruijn Graph Sequence Mapping Problem - https://link.springer.com/chapter/10.1007/978-3-031-36805-9_11
