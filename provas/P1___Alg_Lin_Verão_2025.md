# Primeira Prova - Álgebra Linear (Verão 2025)

A prova é individual e sem consulta, com duração de duas horas. Os resultados dados em sala de aula podem (e devem) ser usados sem demonstração.

---

### Questão 1.
Faça o que se pede.

(a) Discuta a solução do seguinte sistema linear:
```
S = {
  3x + 5y + 12z - w = -3,
  x + y + 4z - w = -6,
  2y + 2z + w = 5.
}
```

(b) Acrescente ao sistema `S` a equação `2z + kw = 9` e encontre um valor de `k` que torne o sistema impossível.

### Questão 2.
Para os dois subespaços vetoriais de R³ que seguem, obtenha um conjunto de vetores geradores.

(a) `U = { (x, y, z) ∈ R³ | x - 2y = 0 }`
(b) `V = { (x, y, z) ∈ R³ | x + 2y - 3z = 0 }`.
(c) Em seguida, calcule os subespaços `U ∩ V` e `U + V`, obtendo bases para os mesmos. Esta última soma é direta?

### Questão 3.
Considere o espaço vetorial sobre R das matrizes quadradas de tamanho 2, `V = M(2, R)`. Sabemos que
`B = { [[1, 0], [0, 0]], [[0, 1], [0, 0]], [[0, 0], [1, 0]], [[0, 0], [0, 1]] }`
é uma possível base (ordenada) para `V`.

(a) Mostre que
`C = { [[1, 0], [0, 2]], [[-1, 0], [0, 0]], [[0, 1], [1, 0]], [[0, 2], [-2, 0]] }`
é também uma base (ordenada) para `V`.

(b) Encontre a matriz mudança da base `B` para a base `C`.

(c) Encontre a matriz mudança da base `C` para a base `B`. Você pode fazer isso tanto pela definição, quanto invertendo a matriz obtida em (b).

(d) Se consideramos a matriz
`A = [[1, 2], [3, 4]] ∈ M(2, R)`,
está claro que na base `B` esta tem coordenadas `(A)B = [[1], [2], [3], [4]]`. Obtenha as coordenadas de `A` na base `C`, isto é `(A)C`.

### Questão 4.
Seja `F : R² → R²` a transformação linear tal que `F(1, 0) = (2, 1)` e `F(0, 1) = (1, 1)`.

(a) Determine `F(2, 4)`.
(b) Determine o vetor `(x, y) ∈ R²` tal que `F(x, y) = (2, 3)`.
(c) Obtenha `ker(F)`, o núcleo de `F`. Este mapa é um isomorfismo (transformação linear bijetora)?
