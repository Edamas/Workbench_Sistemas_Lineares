# Notas para um Curso de Álgebra Linear

**Heitor R. de Assis**

*8 de janeiro de 2026*

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
- 2.1. Axiomas
- 2.2. Bases e dimensão

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

**Definição 1.1.1.** Um corpo F é um conjunto em que duas operações estão definidas: uma chamada de multiplicação `· : F × F → F` e outra chamada de adição `+ : F × F → F`. Estas operações devem satisfazer:
1. `+` é comutativa: `x + y = y + x` para todos `x, y ∈ F`.
2. `+` é associativa: `x + (y + z) = (x + y) + z` para todos `x, y, z ∈ F`.
3. Existe um único elemento neutro `0 ∈ F`, chamado de zero, tal que `x + 0 = 0` para todo `x ∈ F`.
4. Para cada `x ∈ F`, existe um único elemento `(−x) ∈ F` tal que `x + (−x) = 0`.
5. `·` é comutativa: `x · y = y · x` para todos `x, y ∈ F`.
6. `·` é associativa: `x · (y · z) = (x · y) · z` para todos `x, y, z ∈ F`.
7. Existe um único elemento diferente de zero e denotado por 1 que satisfaz `x · 1 = x` para todo `x ∈ F`.
8. Para todo `x ∈ F` diferente de zero, existe um único elemento, denotado por `x⁻¹`, que satisfaz `x⁻¹ · x = 1`.
9. Vale a propriedade distributiva, ou seja, `x · (y + z) = x · y + x · z` para todo `x, y, z ∈ F`.

... (conteúdo completo)
