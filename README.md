# Trabalho de Doutorado

# PMSB — o Problema do Mapeamento de Sequências em Grafo de De Bruijn.

## PMSB

Dadas uma sequência `s` e um grafo de De Bruijn `G`, encontrar o menor percurso `p`
em `G` tal que a sequência induzida por `p` seja a mais parecida possível com `s`.

A tese (`Tese_doutorado.pdf`) trata três variantes desse problema, dependendo de
**onde** as mudanças são permitidas para aproximar `s` de `G`:

| Pasta         | Variante                                             | Onde muda              | Status                                  |
|---------------|-------------------------------------------------------|-------------------------|------------------------------------------|
| `variante1/`  | mudanças na **sequência** (Cap. 5)                    | `s`                     | implementada: algoritmo exato + 4 heurísticas |
| `variante2/`  | mudanças no **grafo** (Cap. 6)                        | `G`                     | implementada: algoritmo exato polinomial (BMTC) |
| `variante3/`  | mudanças em **ambos**, sequência e grafo              | `s` e `G`               | não implementada — problema NP-completo, sem software próprio |

Entre em cada pasta para o `README.md` com detalhes de execução, formato de
entrada/saída e exemplos.

## Painel web (compilar e executar pelo navegador)

```
python3 webui/server.py
```

Abre em `http://localhost:8765` — escolha a ferramenta, compile com um
clique, suba o grafo e a sequência, aperte "Executar" e veja o resultado
(custo destacado quando a ferramenta reporta um). Cada execução fica salva
em `<pasta-da-ferramenta>/results/`. Detalhes e limitações em
`webui/README.md`.

# Publicações
* Heuristics for the de Bruijn Graph Sequence Mapping Problem - https://link.springer.com/chapter/10.1007/978-3-031-36805-9_11

## Autores/Colaboradores

* Lucas Barbosa Rocha
* Francisco Eloi Soares de Araujo
* Said Sadique Adi

## Agradecimentos
* Universidade Federal de Mato Grosso do Sul
* CAPES
