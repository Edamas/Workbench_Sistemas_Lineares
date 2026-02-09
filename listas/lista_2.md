# Segunda Lista de Exercícios - Álgebra Linear (verão 2026)

## Espaços Vetoriais e Transformações Lineares

**Escolhidos para entrega:** P.4 (b) (c), P.5 (a) (b), P.9, P.12, P.14, P.21, P.23 (c) (d), P.27, P.34.

---

### 1. Exercícios Práticos

**P.1** Seja `V = R² = {(x, y) | x, y ∈ R}`. Podemos definir de diversas maneiras as operações soma `+ : R² × R² → R²` e multiplicação por escalar `· : R × R² → R²`. Diga por que as operações definidas em cada item não fazem de `R²` um espaço vetorial.

(a) `(x₁, y₁) + (x₂, y₂) = (x₁ + x₂, 0)`, `λ(x, y) = (λx, λy)`.
(b) `(x₁, y₁) + (x₂, y₂) = (x₁ + x₂, y₁ + y₂)` , `λ(x, y) = (λx, 0)`.
(c) `(x₁, y₁) + (x₂, y₂) = (x₁ + x₂, y₁ + y₂)` , `λ(x, y) = (x, λy)`.
(d) `(x₁, y₁) + (x₂, y₂) = (x₁, y₁)` , `λ(x, y) = (λx, λy)`.
(e) `(x₁, y₁) + (x₂, y₂) = (2x₁ − 2y₁, −x₁ + y₁)` , `λ(x, y) = (3λx, −λx)`.

**P.2** Considere o conjunto `R^∞ = { (x₁, x₂, . . .) | xₖ ∈ R , k ∈ N}`. Mostre que `R^∞`, com as operações `(x₁, x₂, . . .) + (y₁, y₂, . . .) = (x₁ + y₁, x₂ + y₂, . . .)` e `λ · (x₁, x₂, x₃, . . .) = (λx₁, λx₂, λx₃, . . .)` é um espaço vetorial sobre `R`.

**P.3** Sabemos que `V = {f : R → R}` com as operações `(f + g)(x) = f(x) + g(x)` (∀x ∈ R) e `(λf)(x) = λf(x)` (∀x ∈ R) é um espaço vetorial. Agora, mostre que:
(a) o conjunto `C⁰(R, R) = {f : R → R ; f é contínua} ⊂ V` é um subespaço vetorial de V .
(b) generalizando o item (a), o conjunto `Cᵏ(R, R) = {f : R → R ; f possui k-ésima derivada e f⁽ᵏ⁾ é contínua} ⊂ V` é um subespaço vetorial de `V` .
(c) o conjunto `Pₙ(F) = {p : R → R ; p(x) = a₀ + a₁x + · · · + aₙxⁿ}` é um subespaço vetorial de `V`.
(d) Obtenha um subconjunto de `V` que não é subespaço vetorial.

**P.4** Em cada um dos casos, verifique se o conjunto `W` contido no espaço vetorial `V` é um subespaço do espaço vetorial `V`.
(a) `W = {(x, y, z) ∈ F³ | x = 0}`, `V = F³`.
(b) `W = {(x, y, z) ∈ F³ | x ∈ Z}`, `V = F³`.
(c) `W = {p(t) ∈ Pₙ(F) | grau(p) ≥ 2}`, `V = Pₙ(F)`.
(d) `W = {f ∈ C¹(R) | f(t) + f'(t) = 0, ∀t ∈ R}`, `V = C¹(R)`.

... (e assim por diante para todos os exercícios)
