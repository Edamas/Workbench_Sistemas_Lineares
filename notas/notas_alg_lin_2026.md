# Notas para um Curso de Álgebra Linear

**Heitor R. de Assis**

*18 de janeiro de 2026*

---

## Sumário

### 1. Sistemas Lineares e Matrizes
- 1.1. Sistemas Lineares
- 1.2. Resolução de Sistemas Lineares
- 1.3. Matrizes sobre um corpo
- 1.4. Matrizes quadradas
- 1.5. Matrizes e Sistemas Lineares
- 1.6. Algoritmo de inversão de matrizes

### 2. Espaços Vetoriais
- 2.1. Definições
- 2.2. Bases e dimensão
- 2.3. Obtendo bases para subespaços de R^n
- 2.4. Coordenadas e mudança de base

### 3. Transformações Lineares
- 3.1. Transformações lineares e principais propriedades
- 3.2. O espaço das transformações lineares
- 3.3. Teorema do Núcleo e da Imagem
- 3.4. Transformações Lineares e Matrizes

### 4. Espaços Vetoriais com Produto interno
- 4.0.1. Ortogonalidade e ortonormalidade

### 5. Seleção de Exercícios
- 5.1. Sistemas Lineares e Matrizes

---

## Capítulo 1
# Sistemas Lineares e Matrizes

### 1.1 Sistemas Lineares

Sistemas lineares são conjuntos de equações (lineares, portanto o nome) que buscamos resolver simultaneamente. Sua origem é antiga e suas aplicações são as mais variadas, estando presentes em diversos campos da ciência. É a partir do seu estudo que a álgebra linear se formou e, reciprocamente, os sistemas lineares formam uma ferramente extremamente útil para esta teoria.

Os coeficientes presentes nas equações são elementos de um **corpo**. Até o ponto presente, o mais comum é trabalharmos com equações contendo *coeficientes reais*. Se faz importante, porém, equações com coeficientes no corpo dos complexos, C. Pode-se pensar que perguntas sobre os números complexos são tem interesse abstrato, mas a variedade de instâncias onde encontramos números complexos já deve ser suficiente para considerá-los.

**Definição 1.1.1.** Um *corpo* F é um conjunto em que duas operações estão definidas: uma chamada de multiplicação `· : F × F → F` e outra chamada de adição `+ : F × F → F`. Estas operações devem satisfazer:

1. `+` é comutativa: `x + y = y + x` para todos `x, y ∈ F`.
2. `+` é associativa: `x + (y + z) = (x + y) + z` para todos `x, y, z ∈ F`.
3. Existe um único elemento neutro `0 ∈ F`, chamado de zero, tal que `x + 0 = x` para todo `x ∈ F`.
4. Para cada `x ∈ F`, existe um único elemento `(−x) ∈ F` tal que `x + (−x) = 0`.
5. `·` é comutativa: `x · y = y · x` para todos `x, y ∈ F`.
6. `·` é associativa: `x · (y · z) = (x · y) · z` para todos `x, y, z ∈ F`.
7. Existe um único elemento diferente de zero e denotado por 1 que satisfaz `x · 1 = x` para todo `x ∈ F`.
8. Para todo `x ∈ F` diferente de zero, existe um único elemento, denotado por `x⁻¹`, que satisfaz `x⁻¹ · x = 1`.
9. Vale a propriedade distributiva, ou seja, `x · (y + z) = x · y + x · z` para todo `x, y, z ∈ F`.

**Observação.**
Discorramos um pouco mais sobre esta definição, para se ter uma noção de sua significância. Raramente na matemática estuda-se um conjunto de elementos por si só, sendo muito mais comum o estudo de um conjunto de elementos juntamente com as operações que os relacionam. Como exemplo, considere os números reais. Quando falamos da reta real R, não estamos apenas falando do conjunto de números que formam a reta real, `R = {1, 2, 3, π, −2, √3, · · · }`, mas também assumimos que 1+2 = 3, 4π+2π = 6π e assim por diante. Da mesma forma, estamos convencionando que 1×x = x, qualquer valor que `x ∈ R` possa ter, que `−1 × 4 = −4` e que daí `(−1 × 4) + 4 = 0`. Infinitas outras operações entre elementos individuais de R estão de certa forma implícitas quando dizemos estudamos a reta real.

Acontece que é um objetivo do campo da matemática a abstração de conceitos conhecidos. Ao nos depararmos com um conjunto, munido de certas operações que julgamos interessantes, buscamos abstrair o máximo possível. Nos perguntamos, por exemplo: Qual é o conjunto minimal de propriedades que, por si são, são capazes de gerar todas outras? Quando um conjunto, qualquer que seja, satisfaz tais propriedades, quais afirmações podem ser feitas? Por exemplo, a afirmação de que qualquer número multiplicado por zero da zero não está entre as propriedades da Definição 1.1.1. Apesar disso, qualquer conjunto com estas propriedades satisfaz esta afirmação (consegue prová-la?). Esta é uma forma de diferenciarmos as proposições que são intrínsecas de uma certa estrutura daquelas que são específicas de cada modelo. A Definição 1.1.1, portanto, é o resultado da abstração das propriedades essenciais dos números reais.

Esta investigação é útil, inclusive, para testar se as mesmas operações podem ser assumidas sobre diferentes conjuntos. Como exemplo, considere o plano `R² = {(x, y) | x, y ∈ R}`. Temos uma óbvia e bem estabelecida adição neste conjunto, a saber `(x₁, y₁) + (x₂, y₂) = (x₁ + x₂, y₁ + y₂)` . Podemos generalizar a operação de multiplicação que está presente em R para este novo conjunto? A primeira tentativa pode ser replicar o que foi feito para a adição, definindo `(x₁, y₁) × (x₂, y₂) = (x₁ × x₂, y₁ × y₂)` (1.1). Existe alguma propriedade que valorizamos quando falamos de produtos, que R com a multiplicação usual possui, mas R² com a multiplicação definida acima não satisfaz? A resposta é **sim**. Em R, não existem dois números *x* e *y* *diferentes de zero* tais que `x × y = 0`. Já em R², é fácil ver que como definido acima, a multiplicação de duas duplas pode ser zero sem que qualquer uma delas seja zero. Por exemplo, `(x, 0) × (0, y) = (0, 0)`. Acontece que esta é uma propriedade importante por alguns motivos diferentes e, portanto, a definição de produto como em (1.1) não nos é adequada. Podemos nos perguntar se existe alguma definição para o produto em R² que funcione da mesma forma que a multiplicação em R.

Os exemplos mais importantes de corpos são os já conhecidos R e C. Estes não são os únicos, mas serão os únicos considerados aqui. Neste curso, portanto, F será sempre R ou C e sua utilização será principalmente a de um simplificador de notação. Na próxima definição, por exemplo, generalizamos o conceito de produto cartesiano de reais, considerando R e C ao mesmo tempo.

**Definição 1.1.2.** Seja `n ∈ N` um natural qualquer. O conjunto `Fⁿ` é aquele formado por todas as n-uplas de elementos do corpo F, ou seja, `Fⁿ = {(α₁, · · · , αₙ) | α₁, · · · , αₙ ∈ F}`. Neste texto, as letras gregas `α, β, γ`, etc. serão utilizadas para representar elementos do corpo F, seja ele R ou C.

**Definição 1.1.3.** Um *sistema linear* é um conjunto de equações da forma
`(S) { α₁₁x₁ + · · · + α₁ₙxₙ = β₁, ... , αₘ₁x₁ + · · · + αₘₙxₙ = βₘ }`
em que `αᵢⱼ ∈ F` e `β₁, . . . , βₘ ∈ F` são constantes e `x₁, . . . , xₙ` são as incógnitas, cujos valores queremos determinar.

Dizemos que o sistema linear (S) é **homogêneo** se `β₁ = · · · = βₘ = 0`. O conjunto dos elementos `(x₁, . . . , xₙ) ∈ Fⁿ` que resolvem o sistema é chamado de **conjunto solução** do sistema.

Quanto à possibilidade de solução, temos somente três casos:
- Um sistema (S) é **compatível determinado** se existe uma única solução.
- Um sistema (S) é **compatível indeterminado** se existem infinitas soluções.
- Um sistema (S) é **incompatível** se não possui solução.

... (O restante do documento segue a mesma lógica de conversão)
