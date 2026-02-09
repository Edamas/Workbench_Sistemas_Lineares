# GABARITO DA PROVA 1

**Exercício 1.**

Para estudar o sistema, devemos inicialmente escaloná-lo.

```
x + y + z = 0
x + ay + z = 1
x + y - az = a
```

~ L2 - L1

```
x + y + z = 0
(a - 1)y = 1
x + y - az = a
```

~ L3 - L1

```
x + y + z = 0
(a - 1)y = 1
-(a + 1)z = a
```
.

Assim obtemos o sistema escalonado. Ele será incompatível (não terá solução) se, e somente se, obtivermos uma igualdade do tipo `0x + 0y + 0z = β`, em que `β` é diferente de zero. Na expressão acima isto só pode ocorrer nas duas últimas linhas do sistema.

Na segunda linha do sistema escalonado, se `a = 1`, então obtemos uma equação da forma `0y = 1`. Logo o sistema é incompatível, ou seja, não tem solução se `a = 1`.

Na terceira linha do sistema escalonado, se `a = -1`, então obtemos uma equação da forma `0z = -1`. Logo o sistema é incompatível também para `a = -1`.

Se `a` é diferente de `1` e de `-1`, então o sistema escalonado tem 3 equações não nulas e 3 incógnitas. Logo é compatível determinado (tem uma única solução). Para determinar as soluções, podemos proceder de duas maneiras.

**Método 1.**

Da última equação, calculamos que `z = -a / (1+a)`. Da segunda equação achamos que `y = 1 / (a-1)`. Por fim, da primeira achamos que `x = -y - z = 1/(1-a) + a/(1+a) = (1+2a-a²)/(1-a²)`.

**Método 2.**

Continuamos o escalonamento:
```
x + y + z = 0
(a - 1)y = 1
-(a + 1)z = a
```
~ (1/(a-1))L2
~ (-1/(a+1))L3
```
x + y + z = 0
y = 1/(a-1)
z = -a/(a+1)
```
~ L1 - L2
~ L1 - L3
```
x = (1+2a-a²)/(1-a²)
y = 1/(a-1)
z = -a/(-a-1)
```
.

Concluímos que `(x, y, z) = ((1+2a-a²)/(1-a²), 1/(a-1), -a/(1+a))`.

**Exercício 2.**

a) Escrevemos
```
1 = 1 * 1
1 + t = 1 * 1 + 1 * t
1 + t² = 1 * 1 + 1 * t²
t³ = 1 * t³
```
Logo
`ICB = [[1, 1, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]`.

b) Para calcular `IBC`, usamos que `IBC = ICB⁻¹`. Vamos inverter a matriz `ICB`.
```
[ 1 1 1 0 | 1 0 0 0 ]
[ 0 1 0 0 | 0 1 0 0 ]
[ 0 0 1 0 | 0 0 1 0 ]
[ 0 0 0 1 | 0 0 0 1 ]
```
~ L1 - L2
~ L1 - L3
```
[ 1 0 0 0 | 1 -1 -1 0 ]
[ 0 1 0 0 | 0  1  0 0 ]
[ 0 0 1 0 | 0  0  1 0 ]
[ 0 0 0 1 | 0  0  0 1 ]
```
Logo `IBC = [[1, -1, -1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]`.

c) Basta usar `YC = IBC * XB`, em que `XB = (1, 1, 1, 1)`. Neste caso
`[[1, -1, -1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]] * [[1], [1], [1], [1]] = [[-1], [1], [1], [1]]`.
Logo as coordenadas na base C são `(-1, 1, 1, 1)`. Isto pode ser visto também diretamente, pois o polinômio com coordenadas `(1, 1, 1, 1)` na base B é `p(t) = 1 + t + t² + t³`. Este polinômio pode ser escrito como `p(t) = (-1)*1 + 1*(1 + t) + 1*(1 + t²) + 1*t³`.

d) Sim, ambos são isomorfos. Como vimos em sala de aula `M₂(F)` é um espaço vetorial de dimensão 4. `P₃(F)` também é um espaço vetorial de dimensão 4, pois suas bases têm 4 elementos. Usamos agora o resultado que dois espaços vetoriais sobre o mesmo corpo F são isomorfos se, e somente se, têm a mesma dimensão. Como este é o caso, então eles são isomorfos.
Poderíamos resolver este problema também diretamente. De fato, basta definir explicitamente um isomorfismo entre os espaços, tal como `F : P₃(F) → M₂(F)` dado por `F(a + bt + ct² + dt³) = [[a, b], [c, d]]`.
Agora basta provar que a transformação linear dada acima é bijetora, ou seja, injetora e sobrejetora.
É injetora, pois se `F(a + bt + ct² + dt³) = [[0, 0], [0, 0]]`, então `a = b = c = d = 0` e F é injetora.
É sobrejetora, pois seja `[[x, y], [w, z]]` um elemento qualquer de `M₂(F)`. Logo existe um polinômio `x + yt + wt² + zt³`, tal que `F(x + yt + wt² + zt³) = [[x, y], [w, z]]`. Logo a imagem de F é todo o espaço `M₂(F)`, e, portanto, a função é sobrejetora.

**Exercício 3.**

a) A base canônica de `P₂(F)` tem 3 elementos, `1, t` e `t²`. Logo a dimensão, que é o número de elementos que as bases de um certo espaço vetorial finitamente gerado contêm, é igual a 3.

b) Para verificar que `{1 + t², t + t², 1 + t + t²}` é uma base de `P₂(F)`, devemos verificar que o conjunto é linearmente independente (L.I.) e gera o espaço dos polinômios de grau menor ou igual a 2.
Primeiro fato a ser verificado: `{1 + t², t + t², 1 + t + t²}` é um conjunto L.I.
De fato suponha que `α₁(1 + t²) + α₂(t + t²) + α₃(1 + t + t²) = 0`. Rearranjando os termos, obtemos `(α₁ + α₃) + (α₂ + α₃)t + (α₁ + α₂ + α₃)t² = 0`. Como `{1, t, t²}` é base, e, portanto, é linearmente independente, isto só ocorre se `α₁ + α₃ = 0`, `α₂ + α₃ = 0` e `α₁ + α₂ + α₃ = 0`. Ou seja, `α₁(1 + t²) + α₂(t + t²) + α₃(1 + t + t²) = 0` se, e somente se, os `αᵢ` forem solução do sistema linear homogêneo abaixo:
```
α₁ + α₃ = 0
α₂ + α₃ = 0
α₁ + α₂ + α₃ = 0
```
~ L3 - L1
~ L3 - L2
```
α₁ + α₃ = 0
α₂ + α₃ = 0
-α₃ = 0
```
Como o sistema escalonado tem 3 equações não nulas e 3 incógnitas, concluímos que o sistema é compatível determinado, ou seja, tem uma única solução. Como o sistema é homogêneo, concluímos que a única solução é `(α₁, α₂, α₃) = (0, 0, 0)`.
...

**Exercício 4.**
...
