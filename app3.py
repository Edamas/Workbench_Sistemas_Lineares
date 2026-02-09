import streamlit as st
import sympy

# --- FUNÇÕES AUXILIARES DE FORMATAÇÃO E LIMPEZA ---

def cleanup_expr(expr):
    """Remove pequenos artefatos numéricos de uma expressão sympy."""
    if isinstance(expr, sympy.Matrix):
        return expr.applyfunc(cleanup_expr)
    if hasattr(expr, 'evalf'):
        evalf_expr = expr.evalf(20)
        if hasattr(evalf_expr, 'as_real_imag'):
            real_part, imag_part = evalf_expr.as_real_imag()
            if abs(imag_part) < 1e-15:
                real_expr = sympy.re(expr)
                if abs(real_expr.evalf() - round(real_expr.evalf())) < 1e-15:
                    return sympy.Integer(round(real_expr.evalf()))
                return real_expr
    return expr

def format_sympy(expr, notation_style="Fração/Exato"):
    """Formata uma expressão sympy para LaTeX com base no estilo de notação."""
    if expr is None:
        return ""
    
    use_decimal = notation_style == "Decimal"

    if isinstance(expr, sympy.Matrix):
        # Always use pmatrix for matrices
        if use_decimal:
            return sympy.latex(expr.evalf(4), mat_str='pmatrix')
        return sympy.latex(expr, mat_str='pmatrix')
    
    # For non-matrix expressions
    if use_decimal:
        try:
            return sympy.latex(expr.evalf(4))
        except AttributeError:
            pass # Fallback to default if evalf fails
    
    return sympy.latex(expr)

def format_vector_as_tuple_latex(vec, notation_style):
    if not isinstance(vec, sympy.Matrix) or vec.cols != 1:
        return f"\\text{{(Erro: Não é vetor coluna)}}"
    elements_latex = [format_sympy(element, notation_style) for element in vec]
    return f"({', '.join(elements_latex)})"

def format_matrix_as_python_tuple(matrix):
    rows_str = []
    for row_idx in range(matrix.rows):
        row_elements = [str(matrix[row_idx, col_idx]) for col_idx in range(matrix.cols)]
        rows_str.append(f"({', '.join(row_elements)})")
    return f"({', '.join(rows_str)})"
    
# --- FUNÇÃO PRINCIPAL DA APLICAÇÃO ---

def main():
    st.set_page_config(layout="wide", page_title="Calculadora de Diagonalização")
    
    with st.sidebar:
        st.title("Configurações")
        st.markdown("---")
        exemplo_str = "-9 4 4\n-8 3 4\n-16 8 7"
        matrix_str = st.text_area("Insira a matriz A:", exemplo_str, height=100)
        
        notation_style = st.radio("Modo de Notação:", ["Fração/Exato", "Decimal"])
        
        determinant_method = st.radio(
            "Método para det(3x3):",
            ["Expansão de Cofatores", "Regra de Sarrus"]
        )
        
        st.markdown("---")
        run_button = st.button("Calcular Diagonalização", type="primary")

    st.title("Calculadora de Diagonalização de Matrizes")

    if not run_button:
        st.info("Insira a matriz na barra lateral e clique em 'Calcular' para começar.")
        return

    try:
        rows = matrix_str.strip().split('\n')
        matrix_data = [[sympy.Rational(item) for item in row.split()] for row in rows]
        A = sympy.Matrix(matrix_data)

        if not A.is_square:
            st.error("A matriz inserida não é quadrada.")
            return

        st.subheader("Etapa 1: Análise da Matriz de Entrada")
        st.latex(f"A = {format_sympy(A, notation_style)}")

        lambda_ = sympy.Symbol('λ')
        I = sympy.eye(A.rows)
        lambda_I = I * lambda_

        st.subheader("Etapa 2: Encontrando os Autovalores (λ)")

        st.markdown("#### 2.1. Montando o Polinômio Característico, p(λ)")
        det_str = f"p(λ) = \\det \\left( {format_sympy(A, notation_style)} - λ {format_sympy(I, notation_style)} \\right) = \\det \\left( {format_sympy(A, notation_style)} - {format_sympy(lambda_I, notation_style)} \\right)"
        st.latex(det_str)
        
        char_poly_matrix = A - lambda_I
        st.markdown("Calculando a matriz `A - λI`:")
        st.latex(f"A - λI = {format_sympy(char_poly_matrix, notation_style)}")
        
        print(f"DEBUG: Checking A.rows = {A.rows}") # Debug print
        with st.expander("Detalhes do Cálculo do Determinante"):
            if A.rows == 2:
                a, b = char_poly_matrix[0,:]
                c, d = char_poly_matrix[1,:]
                st.markdown("Para uma matriz 2x2, `det = ad - bc`:")
                st.latex(f"p(λ) = ({format_sympy(a, notation_style)})({format_sympy(d, notation_style)}) - ({format_sympy(b, notation_style)})({format_sympy(c, notation_style)})")
            
            elif A.rows == 3 and determinant_method == "Expansão de Cofatores":
                st.markdown("**Usando a Expansão de Cofatores (ao longo da primeira linha):**")
                st.latex(r"\det(A) = a_{11} C_{11} + a_{12} C_{12} + a_{13} C_{13}")
                st.latex(r"C_{ij} = (-1)^{i+j} M_{ij}")
                
                a11 = char_poly_matrix[0,0]
                a12 = char_poly_matrix[0,1]
                a13 = char_poly_matrix[0,2]

                # Minor 11
                M11_matrix = char_poly_matrix.minorMatrix(0,0)
                det_M11 = sympy.det(M11_matrix)
                C11 = ((-1)**(1+1)) * det_M11
                st.markdown("**Cofator $C_{11}$:**")
                st.latex(f"M_{{11}} = {format_sympy(M11_matrix, notation_style)}")
                st.latex(f"\\det(M_{{11}}) = {format_sympy(det_M11, notation_style)}")
                st.latex(f"C_{{11}} = (-1)^{{1+1}} \\cdot \\det(M_{{11}}) = {format_sympy(C11, notation_style)}")

                # Minor 12
                M12_matrix = char_poly_matrix.minorMatrix(0,1)
                det_M12 = sympy.det(M12_matrix)
                C12 = ((-1)**(1+2)) * det_M12
                st.markdown("**Cofator $C_{12}$:**")
                st.latex(f"M_{{12}} = {format_sympy(M12_matrix, notation_style)}")
                st.latex(f"\\det(M_{{12}}) = {format_sympy(det_M12, notation_style)}")
                st.latex(f"C_{{12}} = (-1)^{{1+2}} \\cdot \\det(M_{{12}}) = {format_sympy(C12, notation_style)}")
                
                # Minor 13
                M13_matrix = char_poly_matrix.minorMatrix(0,2)
                det_M13 = sympy.det(M13_matrix)
                C13 = ((-1)**(1+3)) * det_M13
                st.markdown("**Cofator $C_{13}$:**")
                st.latex(f"M_{{13}} = {format_sympy(M13_matrix, notation_style)}")
                st.latex(f"\\det(M_{{13}}) = {format_sympy(det_M13, notation_style)}")
                st.latex(f"C_{{13}} = (-1)^{{1+3}} \\cdot \\det(M_{{13}}) = {format_sympy(C13, notation_style)}")

                # Final Sum
                final_det_cofactor = sympy.expand(a11 * C11 + a12 * C12 + a13 * C13)
                st.markdown("**Soma dos Cofatores Expandidos:**")
                st.latex(f"p(λ) = ({format_sympy(a11, notation_style)})({format_sympy(C11, notation_style)}) + ({format_sympy(a12, notation_style)})({format_sympy(C12, notation_style)}) + ({format_sympy(a13, notation_style)})({format_sympy(C13, notation_style)}) = {format_sympy(final_det_cofactor, notation_style)}")
            
            elif A.rows == 3 and determinant_method == "Regra de Sarrus":
                st.markdown("**Usando a Regra de Sarrus:**")
                sarrus_matrix = char_poly_matrix.col_insert(3, char_poly_matrix.col(0)).col_insert(4, char_poly_matrix.col(1))
                st.latex(f"{format_sympy(sarrus_matrix, notation_style)}")
                
                d1 = char_poly_matrix[0,0] * char_poly_matrix[1,1] * char_poly_matrix[2,2]
                d2 = char_poly_matrix[0,1] * char_poly_matrix[1,2] * char_poly_matrix[2,0]
                d3 = char_poly_matrix[0,2] * char_poly_matrix[1,0] * char_poly_matrix[2,1]
                sum_main = sympy.expand(d1 + d2 + d3)

                a1 = char_poly_matrix[0,2] * char_poly_matrix[1,1] * char_poly_matrix[2,0]
                a2 = char_poly_matrix[0,0] * char_poly_matrix[1,2] * char_poly_matrix[2,1]
                a3 = char_poly_matrix[0,1] * char_poly_matrix[1,0] * char_poly_matrix[2,2]
                sum_anti = sympy.expand(a1 + a2 + a3)

                st.markdown("**Soma dos produtos das diagonais principais:**")
                st.latex(f"({format_sympy(sympy.expand(d1), notation_style)}) + ({format_sympy(sympy.expand(d2), notation_style)}) + ({format_sympy(sympy.expand(d3), notation_style)}) = {format_sympy(sum_main, notation_style)}")
                st.markdown("**Soma dos produtos das diagonais secundárias:**")
                st.latex(f"({format_sympy(sympy.expand(a1), notation_style)}) + ({format_sympy(sympy.expand(a2), notation_style)}) + ({format_sympy(sympy.expand(a3), notation_style)}) = {format_sympy(sum_anti, notation_style)}")
                st.markdown("**Determinante (Soma Principais - Soma Secundárias):**")
                st.latex(f"p(λ) = ({format_sympy(sum_main, notation_style)}) - ({format_sympy(sum_anti, notation_style)})")

        char_poly = A.charpoly(lambda_).as_expr()
        st.markdown("O polinômio característico final é:")
        st.latex(rf"p(λ) = {format_sympy(char_poly, notation_style)} = 0")
        
        st.markdown("#### 2.2. Encontrando as Raízes do Polinômio")
        with st.expander("Detalhes da Resolução do Polinômio"):
            poly = sympy.Poly(char_poly, lambda_)
            
            if poly.degree() == 2:
                st.markdown("Para um polinômio de grau 2, aplicamos a Fórmula de Bhaskara:")
                st.latex(r"λ = \frac{{-b \pm \sqrt{{b^2 - 4ac}}}}{{2a}}")
                a, b, c = poly.all_coeffs()
                st.markdown(f"Para `{format_sympy(poly.as_expr(), notation_style)}`, temos `a={format_sympy(a, notation_style)}`, `b={format_sympy(b, notation_style)}`, `c={format_sympy(c, notation_style)}`.")
                delta = b**2 - 4*a*c
                st.latex(rf"\Delta = b^2 - 4ac = {format_sympy(delta, notation_style)}")
                if delta >= 0:
                    lambda1 = (-b + sympy.sqrt(delta)) / (2*a)
                    lambda2 = (-b - sympy.sqrt(delta)) / (2*a)
                    st.markdown("As raízes são:")
                    st.latex(rf"λ_1 = {format_sympy(lambda1, notation_style)}, \; λ_2 = {format_sympy(lambda2, notation_style)}")
                else:
                    st.markdown("O discriminante é negativo, indicando raízes complexas. As raízes serão apresentadas abaixo.")

            elif poly.degree() == 3:
                st.markdown("**Para um polinômio de grau 3, buscamos por raízes racionais (Teorema das Raízes Racionais).**")
                
                coeffs = poly.all_coeffs()
                leading_coeff = coeffs[0]
                constant_term = coeffs[-1]

                divisors_p = sympy.divisors(abs(constant_term))
                divisors_q = sympy.divisors(abs(leading_coeff))
                
                possible_roots_str = []
                possible_roots_set = set()
                for p_val in divisors_p:
                    for q_val in divisors_q:
                        root1 = sympy.Rational(p_val, q_val)
                        root2 = sympy.Rational(-p_val, q_val)
                        
                        if root1 not in possible_roots_set:
                            possible_roots_str.append(format_sympy(root1, notation_style))
                            possible_roots_set.add(root1)
                        if root2 not in possible_roots_set:
                            possible_roots_str.append(format_sympy(root2, notation_style))
                            possible_roots_set.add(root2)
                
                st.markdown(f"1. **Identificar `p` e `q`:** O termo constante é `{format_sympy(constant_term, notation_style)}` e o coeficiente líder é `{format_sympy(leading_coeff, notation_style)}`.")
                st.markdown(f"2. **Divisores de `p`:** `{', '.join([format_sympy(d, notation_style) for d in divisors_p])}`")
                st.markdown(f"3. **Divisores de `q`:** `{', '.join([format_sympy(d, notation_style) for d in divisors_q])}`")
                st.markdown(f"4. **Possíveis raízes racionais (p/q):** `{', '.join(possible_roots_str)}`")
                st.markdown("5. **Testamos estas raízes no polinômio.**")

                # Implementação manual do teste de raízes racionais
                r = None
                for root_candidate in possible_roots_set:
                    if poly.subs(lambda_, root_candidate) == 0:
                        r = root_candidate
                        break

                if r is not None:
                    st.markdown(f"Um teste revela uma raiz racional em `λ₁ = {format_sympy(r, notation_style)}`.")
                    st.markdown(f"Dividindo o polinômio por `(λ - ({format_sympy(r, notation_style)}))`, obtemos um polinômio quadrático:")
                    quad_poly, remainder = sympy.div(poly, (lambda_ - r))
                    st.latex(f"{format_sympy(poly.as_expr(), notation_style)} / (λ - {format_sympy(r, notation_style)}) = {format_sympy(quad_poly.as_expr(), notation_style)}")
                    
                    st.markdown("Agora, aplicamos a Fórmula de Bhaskara para encontrar as raízes restantes:")
                    st.latex(r"λ = \frac{{-b \pm \sqrt{{b^2 - 4ac}}}}{{2a}}")
                    a, b, c = quad_poly.all_coeffs()
                    st.markdown(f"Para `{format_sympy(quad_poly.as_expr(), notation_style)}`, temos `a={format_sympy(a, notation_style)}`, `b={format_sympy(b, notation_style)}`, `c={format_sympy(c, notation_style)}`.")
                    delta = b**2 - 4*a*c
                    st.latex(rf"\Delta = b^2 - 4ac = {format_sympy(delta, notation_style)}")
                    lambda2 = (-b + sympy.sqrt(delta)) / (2*a)
                    lambda3 = (-b - sympy.sqrt(delta)) / (2*a)
                    st.markdown("As outras raízes são:")
                    st.latex(rf"λ_2 = {format_sympy(lambda2, notation_style)}, \; λ_3 = {format_sympy(lambda3, notation_style)}")
                else:
                    st.markdown("Não foi encontrada uma raiz racional exata por este método de teste.")
                    st.markdown("As raízes são encontradas computacionalmente.")

        eigenvals_dict = sympy.roots(char_poly, lambda_)
        autovalores = sorted(eigenvals_dict.keys(), key=lambda x: x.evalf())
        cleaned_autovalores = [cleanup_expr(v) for v in autovalores]
        
        st.markdown("Os autovalores da matriz A são:")
        st.latex(rf"λ \in {format_sympy(sympy.FiniteSet(*cleaned_autovalores), notation_style)}")
        
        cleaned_eigenvals_dict = {cleanup_expr(k): v for k, v in eigenvals_dict.items()}
        st.markdown("**Análise dos Autovalores:**")
        for val, mult in cleaned_eigenvals_dict.items():
            st.write(f"- O autovalor **{format_sympy(val, notation_style)}** tem multiplicidade algébrica **{mult}**.")
        
        st.subheader("Etapa 3: Encontrando os Autovetores (v)")
        all_eigenvectors_map = {}
        for val in autovalores:
            cleaned_val = cleanup_expr(val)
            basis = (A - cleaned_val * I).nullspace()
            all_eigenvectors_map[cleaned_val] = [cleanup_expr(vec) for vec in basis]
        
        for cleaned_val, cleaned_basis in all_eigenvectors_map.items():
            st.markdown(f"#### 3.1. Autoespaço para λ = {format_sympy(cleaned_val, notation_style)}")
            st.latex(f"\\text{{Resolvendo }}(A - ({format_sympy(cleaned_val, notation_style)})I)v = 0")
            
            with st.expander(f"Detalhes do Cálculo de A - ({format_sympy(cleaned_val, notation_style)})I"):
                lambda_sym = sympy.Symbol('λ') # Define lambda_sym to represent lambda here
                I_val = sympy.eye(A.rows)
                lambda_I_matrix = I_val * cleaned_val
                
                st.markdown(f"Calculando `({format_sympy(cleaned_val, notation_style)})I`:")
                st.latex(f"({format_sympy(cleaned_val, notation_style)})I = {format_sympy(I_val, notation_style)} \\cdot {format_sympy(cleaned_val, notation_style)} = {format_sympy(lambda_I_matrix, notation_style)}")
                
                char_matrix_eigen = A - lambda_I_matrix
                st.markdown(f"Calculando `A - ({format_sympy(cleaned_val, notation_style)})I`:")
                st.latex(f"A - ({format_sympy(cleaned_val, notation_style)})I = {format_sympy(A, notation_style)} - {format_sympy(lambda_I_matrix, notation_style)} = {format_sympy(char_matrix_eigen, notation_style)}")
                
                # Optional: Show row reduction of char_matrix_eigen
                st.markdown("Reduzindo a matriz à forma escalonada:")
                rref_matrix, pivots = char_matrix_eigen.rref()
                st.latex(f"\\text{{rref}}(A - ({format_sympy(cleaned_val, notation_style)})I) = {format_sympy(rref_matrix, notation_style)}")

            basis_latex = ", ".join([format_vector_as_tuple_latex(vec, notation_style) for vec in cleaned_basis])
            st.markdown("A base para o espaço nulo (kernel) é:")
            st.latex(rf"E_λ({format_sympy(cleaned_val, notation_style)}) = \text{{gerado por}} \left\{{ {basis_latex} \right\}}")
        
        is_diagonalizable = all(cleaned_eigenvals_dict.get(v, 0) == len(all_eigenvectors_map.get(v, [])) for v in cleaned_eigenvals_dict)

        st.subheader("Etapa 4: Construindo as Matrizes P e D")
        if not is_diagonalizable:
            st.error("A matriz A não é diagonalizável.")
            return

        st.success("A matriz A é diagonalizável!")
        
        sorted_autovalores = sorted(cleaned_eigenvals_dict.keys(), key=lambda x: x.evalf())
        P_cols, D_vals = [], []
        for val in sorted_autovalores:
            mult = cleaned_eigenvals_dict[val]
            D_vals.extend([val] * mult)
            P_cols.extend(all_eigenvectors_map.get(val, []))

        P, D = sympy.Matrix.hstack(*P_cols), sympy.diag(*D_vals)
        st.markdown("A matriz **P** é formada pelos autovetores:")
        st.latex(f"P = {format_sympy(P, notation_style)}")
        st.write(f"P (formato textual): {format_matrix_as_python_tuple(P)}")
        st.markdown("E **D** pelos autovalores:")
        st.latex(f"D = {format_sympy(D, notation_style)}")
        
        P_inv = P.inv()
        st.markdown("A inversa de P é:")
        st.latex(f"P^{{-1}} = {format_sympy(P_inv, notation_style)}")
        st.write(f"P⁻¹ (formato textual): {format_matrix_as_python_tuple(P_inv)}")

        st.subheader("Etapa 5: Verificação Final")
        result_matrix = P * D * P_inv
        final_result = result_matrix
        
        st.latex(f"A = PDP^{{-1}}")
        st.latex(f"{format_sympy(P, notation_style)}{format_sympy(D, notation_style)}{format_sympy(P_inv, notation_style)}")
        st.latex(f" = {format_sympy(final_result, notation_style)}")

        if (A - final_result).is_zero_matrix:
            st.success("Verificação confirmada: A = PDP⁻¹ é verdade!")
        else:
            st.warning("A verificação A = PDP⁻¹ falhou.")

    except Exception as e:
        st.error(f"Ocorreu um erro durante o cálculo: {e}")

if __name__ == "__main__":
    main()