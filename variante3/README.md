# variante3 — Mapeamento de Sequências em Grafo de De Bruijn (mudanças na sequência E no grafo)

## O problema

Esta é a **variante 3** do PMSB, mencionada na tese (`../Tese_doutorado.pdf`)
junto com as variantes 1 (`../variante1/`, muda só a sequência) e 2
(`../variante2/`, muda só o grafo): aqui, **tanto a sequência `s` quanto o
grafo de De Bruijn `G` podem ser alterados** para que exista um percurso em
`G` que induza (uma versão editada de) `s`, minimizando o custo total das
duas edições combinadas.

## Por que não há software aqui

Essa variante combinada não foi implementada. Seguindo a classificação de
Amir et al. (que a tese usa como referência para as três variantes do
problema de mapeamento aproximado em hipertexto — Cap. 2/3), permitir edição
simultânea do padrão e da estrutura (variante 3, `PMAPHp3`/`PMSBp3`) cai no
caso que os autores **provaram ser NP-completo**, diferente da variante 1
(polinomial, Amir et al.) e da variante 2 tratada nesta tese (também
polinomial, graças à indução de arcos discutida em `../variante2/README.md`).

Como o problema em si é NP-completo nesta formulação — e a tese não propõe
(nem tenta provar) uma solução polinomial pra ela — não existe aqui um
algoritmo exato equivalente ao `bmt`/`btl` das outras variantes, nem faria
sentido implementar uma heurística sem antes ter uma formulação e uma análise
de complexidade próprias para o problema combinado. Esta pasta fica como
registro de que a variante 3 é conhecida e foi conscientemente deixada de
fora do escopo da tese — não é um esquecimento, é um problema em aberto
(ver `todo`).

## Trabalho futuro

Se algum dia isso for retomado, os caminhos naturais são: (a) restringir o
espaço de busca (ex.: limitar o número de edições permitidas em `s` e em
`G` separadamente) para tentar uma solução exata em tempo pseudo-polinomial;
ou (b) desenvolver uma heurística inspirada nas de `../variante1/bmt/`
(seed-and-extend) combinada com a edição de grafo de `../variante2/btl/`,
sem garantia de otimalidade.
