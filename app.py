import streamlit as st
import sympy
from sympy import latex, Matrix
import system_generator
import numpy as np 
import re 

# Import parsing functions from parser.py (the new, robust version)
from parser import parse_system_input, parse_to_sympy

# --- Configuração da Página ---
st.set_page_config(layout="wide", page_title="Workbench de Sistemas Lineares")

# --- Funções de Lógica e Análise ---

def get_pivot_info(matrix):
    pivots = {}
    zero_rows = []
    contradiction_rows = []
    if not hasattr(matrix, 'shape'):
        return {'pivots': pivots, 'zero_rows': zero_rows, 'contradiction_rows': contradiction_rows}
    num_rows, num_cols = matrix.shape
    for r_idx in range(num_rows):
        row = matrix.row(r_idx)
        is_coeffs_zero = all(sympy.sympify(c).is_zero for c in row[:-1])
        if is_coeffs_zero:
            if sympy.sympify(row[-1]).is_zero:
                zero_rows.append(r_idx)
            else:
                contradiction_rows.append(r_idx)
        else:
            for c_idx in range(num_cols - 1):
                if not sympy.sympify(row[c_idx]).is_zero:
                    pivots[r_idx] = c_idx
                    break
    return {'pivots': pivots, 'zero_rows': zero_rows, 'contradiction_rows': contradiction_rows}

def analyze_row_factors(matrix, row_idx, pivot_info):
    factors = { 'tipo': ('➖', 'Indefinido'), 'unitario': ('➖', 'Não aplicável'), 'reduzida': ('➖', 'Não aplicável') }
    if row_idx in pivot_info['pivots']:
        factors['tipo'] = ('✅', 'Linha de Pivô')
    elif row_idx in pivot_info['contradiction_rows']:
        factors['tipo'] = ('❗️', 'Linha de Contradição (0=k)')
    elif row_idx in pivot_info['zero_rows']:
        factors['tipo'] = ('✳️', 'Linha Nula (0=0)')
    if factors['tipo'][0] == '✅':
        pivot_col = pivot_info['pivots'][row_idx]
        pivot_val = matrix[row_idx, pivot_col]
        factors['unitario'] = ('✅', 'Pivô é 1') if pivot_val == 1 else ('⚠️', 'Pivô não é 1')
        num_vars = matrix.shape[1] - 1
        is_row_reduced = all(sympy.sympify(matrix[row_idx, c]).is_zero for c in range(num_vars) if c != pivot_col)
        factors['reduzida'] = ('✅', 'Linha Reduzida') if is_row_reduced else ('❌', 'Linha Não-Reduzida')
    return factors

def analyze_col_factors(matrix, col_idx, pivot_info, num_vars):
    factors = { 'tipo': ('➖', 'Indefinido'), 'unitario': ('➖', 'Não aplicável'), 'reduzida': ('➖', 'Não aplicável') }
    pivot_cols = {v: k for k, v in pivot_info['pivots'].items()}
    if col_idx < num_vars:
        factors['tipo'] = ('✅', 'Coluna de Pivô') if col_idx in pivot_cols else ('ℹ️', 'Coluna Livre')
    else:
        if pivot_info['contradiction_rows']: factors['tipo'] = ('❗️', 'Coluna de Contradição')
        elif pivot_info['zero_rows'] and not pivot_info['pivots']: factors['tipo'] = ('✳️', 'Coluna Nula')
        else: factors['tipo'] = ('➡️', 'Coluna de Resultados')
    if factors['tipo'][0] == '✅':
        pivot_row = pivot_cols[col_idx]
        pivot_val = matrix[pivot_row, col_idx]
        factors['unitario'] = ('✅', 'Pivô é 1') if pivot_val == 1 else ('⚠️', 'Pivô não é 1')
        is_col_reduced = all(sympy.sympify(matrix[r, col_idx]).is_zero for r in range(matrix.shape[0]) if r != pivot_row)
        factors['reduzida'] = ('✅', 'Coluna Reduzida') if is_col_reduced else ('❌', 'Coluna Não-Reduzida')
    return factors

def calculate_rref_progress(matrix, variables):
    num_rows, num_vars = matrix.shape[0], len(variables)
    if num_rows == 0: return 0.0
    pivot_info = get_pivot_info(matrix)
    rref_compliant_rows = 0
    for i in range(num_rows):
        row_analysis = analyze_row_factors(matrix, i, pivot_info)
        row_type = row_analysis['tipo'][0]
        if row_type == '✳️':
            rref_compliant_rows += 1
            continue
        if row_type == '✅':
            is_pivot_one = row_analysis['unitario'][0] == '✅'
            pivot_col_idx = pivot_info['pivots'][i]
            col_analysis = analyze_col_factors(matrix, pivot_col_idx, pivot_info, num_vars)
            is_col_reduced = col_analysis['reduzida'][0] == '✅'
            if is_pivot_one and is_col_reduced:
                rref_compliant_rows += 1
    return rref_compliant_rows / num_rows if num_rows > 0 else 0.0

def classify_system(matrix, num_vars):
    rank_A = matrix[:,:-1].rank()
    rank_Ab = matrix.rank()
    if rank_A != rank_Ab:
        return "Sistema Impossível (SI)", f"posto(A) = {rank_A} ≠ posto(A|b) = {rank_Ab}"
    else:
        if rank_A == num_vars: return "Sistema Possível e Determinado (SPD)", f"posto(A) = posto(A|b) = n = {rank_A}"
        else: return "Sistema Possível e Indeterminado (SPI)", f"posto(A) = posto(A|b) = {rank_A} < n = {num_vars}"

def solve_echelon_spd(echelon_matrix, variables):
    try:
        sympy_vars = [sympy.Symbol(v) for v in variables]
        solution = sympy.solve_linear_system(echelon_matrix, *sympy_vars)
        if solution is None: return None, "O sistema não possui solução única."
        return {str(var): solution.get(var) for var in sympy_vars}, None
    except Exception as e: return None, f"Erro ao resolver: {e}"

def is_row_echelon(matrix):
    pivots = get_pivot_info(matrix)['pivots']
    last_pivot_col = -1
    for r_idx in sorted(pivots.keys()):
        pivot_col = pivots[r_idx]
        if pivot_col <= last_pivot_col: return False
        for below_r_idx in range(r_idx + 1, matrix.shape[0]):
            if not matrix[below_r_idx, pivot_col].is_zero: return False
        last_pivot_col = pivot_col
    return True

def is_reduced_row_echelon(matrix):
    if not is_row_echelon(matrix): return False
    pivots = get_pivot_info(matrix)['pivots']
    for r_idx, pivot_col in pivots.items():
        if not matrix[r_idx, pivot_col] == 1: return False
        for i in range(matrix.shape[0]):
            if i != r_idx and not matrix[i, pivot_col].is_zero: return False
    return True

# --- Funções de Automação ---

def find_next_ref_step(matrix):
    if is_row_echelon(matrix): return None
    
    num_rows, num_cols = matrix.shape
    pivot_row = 0
    for j in range(num_cols - 1):
        if pivot_row >= num_rows: break
        
        i = pivot_row
        if matrix[i, j] == 0:
            for k in range(i + 1, num_rows):
                if matrix[k, j] != 0:
                    return {'op_type': 'Troca de Linhas', 'r1': k + 1, 'r2': i + 1}
        
        if matrix[pivot_row, j] != 0:
            for i in range(pivot_row + 1, num_rows):
                if matrix[i, j] != 0:
                    factor = -matrix[i, j] / matrix[pivot_row, j]
                    return {'op_type': 'Combinação Linear', 'ri': i + 1, 'rj': pivot_row + 1, 'factor_lambda': str(factor)}
            pivot_row += 1
    return None

def find_next_step(matrix):
    if is_reduced_row_echelon(matrix): return None 
    
    ref_step = find_next_ref_step(matrix)
    if ref_step: return ref_step

    pivots = get_pivot_info(matrix)['pivots']
    for r, c in sorted(pivots.items(), reverse=True):
        if matrix[r, c] != 1:
            factor = 1 / matrix[r, c]
            return {'op_type': 'Multiplicação por Escalar', 'r': r + 1, 'alpha': str(factor)}
        for i in range(r):
            if matrix[i, c] != 0:
                factor = -matrix[i, c]
                return {'op_type': 'Combinação Linear', 'ri': i + 1, 'rj': r + 1, 'factor_lambda': str(factor)}
    return None

def apply_next_step():
    op = find_next_step(st.session_state.current_matrix)
    if op:
        op_type = op.pop('op_type')
        apply_operation(op_type, **op, source="auto")
    else:
        st.sidebar.success("A matriz já está na RREF!")

def solve_to_ref():
    MAX_STEPS = 50 
    for i in range(MAX_STEPS):
        op = find_next_ref_step(st.session_state.current_matrix)
        if op:
            op_type = op.pop('op_type')
            apply_operation(op_type, **op, source="auto")
        else:
            st.sidebar.success("Sistema escalonado (REF)!")
            return
    st.sidebar.error(f"A resolução parou após {MAX_STEPS} passos.")

def solve_all_steps():
    MAX_STEPS = 50 
    for _ in range(MAX_STEPS):
        op = find_next_step(st.session_state.current_matrix)
        if op:
            op_type = op.pop('op_type')
            apply_operation(op_type, **op, source="auto")
        else:
            st.sidebar.success("Sistema resolvido para RREF!")
            return
    st.sidebar.error(f"A resolução automática parou após {MAX_STEPS} passos.")

# --- Funções de UI e Formatação ---

def render_badge(icon, tooltip):
    st.markdown(f"<div style='text-align: center;' title='{tooltip}'>{icon}</div>", unsafe_allow_html=True)

def to_subscript(s):
    match = re.match(r'([a-zA-Z]+)(\d+)', s)
    if not match: return s
    name, digits = match.groups()
    subscript_digits = "".join(['₀₁₂₃₄₅₆₇₈₉'[int(d)] for d in digits])
    return name + subscript_digits

def number_to_subscript(n):
    return "".join(['₀₁₂₃₄₅₆₇₈₉'[int(d)] for d in str(n)])

def format_plain_text(expr):
    if expr is None: return ""
    try:
        if expr.is_Integer: return str(expr)
        if expr.is_Float:
            val = float(expr)
            if np.isclose(val, round(val)): return str(int(round(val)))
            return f"{val:.3f}".rstrip('0').rstrip('.')
        if expr.is_Rational: return str(expr)
    except (AttributeError, TypeError): pass
    return sympy.pretty(expr, use_unicode=True)

def to_latex_equation(row_vector, variables, is_absurd=False):
    line = []
    for j, var_name in enumerate(variables):
        coeff = row_vector[j]
        if not sympy.sympify(coeff).is_zero:
            term_latex = latex(coeff * sympy.Symbol(var_name), mul_symbol=' ')
            if coeff == 1: term_latex = latex(sympy.Symbol(var_name))
            elif coeff == -1: term_latex = f"-{latex(sympy.Symbol(var_name))}"
            line.append(term_latex if not line else (f"+ {term_latex}" if not term_latex.startswith('-') else term_latex))
    lhs = " ".join(line).replace("+ -", "- ").strip() or "0"
    op = r"\\neq" if is_absurd else "="
    return f"{lhs} {op} {latex(row_vector[-1])}"

def to_latex_system(matrix, variables):
    begin_system = r"S = \begin{cases}"
    end_system = r"\end{cases}"
    content = r" \\ ".join([to_latex_equation(list(matrix.row(i)), variables) for i in range(matrix.shape[0])])
    return begin_system + content + end_system

def to_plain_text_equation(row_vector, variables, is_absurd=False):
    line = []
    for j, var_name in enumerate(variables):
        coeff = row_vector[j]
        if not sympy.sympify(coeff).is_zero:
            term_text = f"{format_plain_text(coeff)}*{to_subscript(var_name)}"
            if coeff == 1: term_text = to_subscript(var_name)
            elif coeff == -1: term_text = f"-{to_subscript(var_name)}"
            if not line or term_text.startswith("-"): line.append(f" {term_text}")
            else: line.append(f" + {term_text}")
    lhs = "".join(line).strip() or "0"
    op = " ≠ " if is_absurd else " = "
    return f"{lhs}{op}{format_plain_text(row_vector[-1])}"

def to_plain_text_system(matrix, variables):
    return "\n".join([to_plain_text_equation(list(matrix.row(i)), variables) for i in range(matrix.shape[0])])

def generate_solution_log_items():
    if not st.session_state.get('matrix_history'): return []
    
    def compact_format(matrix):
        formatted_rows = []
        for row_idx in range(matrix.rows):
            formatted_elements = [format_plain_text(elem) for elem in matrix.row(row_idx)]
            formatted_rows.append(f"[{', '.join(formatted_elements)}]")
        return "\n".join(formatted_rows)
    
    items = []
    original_input_str = st.session_state.get('user_input', '')
    if original_input_str:
        items.append(('header', "**Sistema de Equações Original:**"))
        items.append(('matrix_str', original_input_str))

    if st.session_state.get('matrix_history'):
        items.append(('header', "**Matriz Aumentada Inicial [A|b]:**"))
        items.append(('matrix_str', compact_format(st.session_state.matrix_history[0])))

    if not st.session_state.history:
        initial_matrix = st.session_state.matrix_history[0]
        if is_reduced_row_echelon(initial_matrix):
            items.append(('header', "**Observação:** Nenhuma operação foi necessária."))
    else:
        for i, (plain_op, _) in enumerate(st.session_state.history):
            if plain_op: 
                items.append(('header', f"**Passo {i+1}:** `{plain_op}`"))
                items.append(('matrix_str', compact_format(st.session_state.matrix_history[i+1])))
    
    items.append(('header', "---\n\n**ANÁLISE FINAL DO SISTEMA**\n"))
    final_matrix = st.session_state.current_matrix
    variables = st.session_state.variables
    num_vars = len(variables)
    
    final_matrix_label = "**1. Matriz Final:**"
    items.append(('header', final_matrix_label))
    items.append(('matrix_str', compact_format(final_matrix)))

    rank_A, rank_Ab = final_matrix[:,:-1].rank(), final_matrix.rank()
    analysis_text = f"**2. Análise dos Postos:**\n- Posto(A) = {rank_A}\n- Posto(A|b) = {rank_Ab}\n- N° de incógnitas = {num_vars}\n\n**3. Conclusão (Teorema de Rouché-Capelli):**"
    items.append(('header', analysis_text))
    classification, _ = classify_system(final_matrix, num_vars)

    if classification == "Sistema Impossível (SI)":
        items.append(('header', f"Como posto(A) ≠ posto(A|b), o sistema é **Impossível (SI)**.\n- **Conjunto Solução (S):** S = ∅"))
    elif classification == "Sistema Possível e Determinado (SPD)":
        items.append(('header', f"Como posto(A) = posto(A|b) = n, o sistema é **Possível e Determinado (SPD)**."))
        solution, error = solve_echelon_spd(final_matrix, variables)
        if not error:
            solution_text = ", ".join([f"{to_subscript(var)} = {format_plain_text(val)}" for var, val in solution.items() if val is not None])
            if solution_text:
                items.append(('header', f"\n**4. Obtenção da Solução:**\n- **Conjunto Solução (S):** S = {{{solution_text}}}"))
    else:
        graus_liberdade = num_vars - rank_A
        items.append(('header', f"Como posto(A) = posto(A|b) = {rank_A} < n, o sistema é **Possível e Indeterminado (SPI)** com **{graus_liberdade}** grau(s) de liberdade."))
        items.append(('header', "\n**4. Obtenção da Solução Geral:**"))
        items.append(('matrix_str', to_plain_text_system(final_matrix, variables)))
    return items

# --- Lógica de Sessão e Operações Manuais ---
if 'system_processed' not in st.session_state: st.session_state.system_processed = False
if 'history' not in st.session_state: st.session_state.update({'history':[], 'matrix_history':[], 'redo_history':[], 'redo_matrix_history':[]})
if 'run_id' not in st.session_state: st.session_state.run_id = 0
if 'user_input' not in st.session_state: st.session_state.user_input = "2x+y-z=8\n-3x-y+2z=-11\n-2x+y+2z=-3"

def apply_operation(op_type, source="manual", **kwargs):
    if source == "manual":
        st.session_state.redo_history = []
        st.session_state.redo_matrix_history = []
    
    matrix = st.session_state.current_matrix.copy()
    try:
        plain_desc, latex_desc = "", ""
        if op_type == "Troca de Linhas":
            r1, r2 = kwargs['r1'] - 1, kwargs['r2'] - 1
            matrix.row_swap(r1, r2)
            plain_desc = f"L{r1+1} <-> L{r2+1}"
            latex_desc = rf"L_{{{r1+1}}} \leftrightarrow L_{{{r2+1}}}"
        elif op_type == "Multiplicação por Escalar":
            r, alpha_str = kwargs['r'] - 1, kwargs['alpha']
            alpha = parse_to_sympy(alpha_str)
            if alpha == 0: raise ValueError("O escalar (α) não pode ser zero.")
            for j in range(matrix.cols): matrix[r, j] = sympy.simplify(matrix[r, j] * alpha)
            plain_desc = f"L{r+1} <- ({alpha_str}) * L{r+1}"
            latex_desc = rf"L_{{{r+1}}} \leftarrow ({latex(alpha)}) \cdot L_{{{r+1}}}"
        elif op_type == "Combinação Linear":
            ri, rj, lambda_str = kwargs['ri'] - 1, kwargs['rj'] - 1, kwargs['factor_lambda']
            lambda_val = parse_to_sympy(lambda_str)
            if lambda_val == 0: raise ValueError("O fator (λ) não pode ser zero.")
            for j in range(matrix.cols): matrix[ri, j] = sympy.simplify(matrix[ri, j] + lambda_val * matrix[rj, j])
            plain_desc = f"L{ri+1} <- L{ri+1} + ({lambda_str}) * L{rj+1}"
            latex_desc = rf"L_{{{ri+1}}} \leftarrow L_{{{ri+1}}} + ({latex(lambda_val)}) \cdot L_{{{rj+1}}}"
        else: return

        op_type_abbr = {"Troca de Linhas": "Troca", "Multiplicação por Escalar": "Mult.", "Combinação Linear": "Comb. Lin."}.get(op_type, op_type)
        plain_desc += f" [{op_type_abbr}]"

        if is_reduced_row_echelon(matrix):
            plain_desc += " [Matriz na Forma Escalonada Reduzida]"
        elif is_row_echelon(matrix):
            plain_desc += " [Matriz na Forma Escalonada]"

        st.session_state.current_matrix = matrix
        st.session_state.history.append((plain_desc, latex_desc))
        st.session_state.matrix_history.append(matrix)

        if source == "manual":
            st.sidebar.success("Operação aplicada!")
            st.sidebar.latex(latex_desc)

    except Exception as e:
        st.sidebar.error(f"Erro na operação: {e}")
        st.session_state.run_id += 1

def undo_last_operation():
    if len(st.session_state.history) > 0:
        st.session_state.redo_history.append(st.session_state.history.pop())
        st.session_state.redo_matrix_history.append(st.session_state.matrix_history.pop())
        st.session_state.current_matrix = st.session_state.matrix_history[-1]

def redo_next_operation():
    if len(st.session_state.redo_history) > 0:
        st.session_state.history.append(st.session_state.redo_history.pop())
        st.session_state.matrix_history.append(st.session_state.redo_matrix_history.pop())
        st.session_state.current_matrix = st.session_state.matrix_history[-1]

# --- Interface Streamlit (Frontend) ---

with st.sidebar:
    st.title("Workbench de Sistemas Lineares")
    
    with st.expander("0. Gerador de Sistema", expanded=False):
        # Using new keys to avoid conflicts if main app uses them
        gen_num_vars = st.slider("Número de variáveis", 1, 10, 3, key='gen_num_vars_sidebar')
        gen_var_notation = st.radio("Notação das variáveis", ["x, y, z...", "x₁, x₂, x₃..."], index=1, key='gen_var_notation_sidebar')
        
        c1, c2, c3 = st.columns(3)
        gen_use_fractions = c1.checkbox("Frações", value=False, key='gen_use_fractions_sidebar')
        gen_use_complex = c2.checkbox("Complexos", key='gen_use_complex_sidebar')
        gen_use_roots = c3.checkbox("Raízes", key='gen_use_roots_sidebar')
        
        gen_solution_type = st.selectbox("Tipo de Solução", ["Determinado", "Indeterminado", "Impossível"], key='gen_solution_type_sidebar')
        is_homogeneous = st.checkbox("Sistema homogêneo", value=False, key='is_homogeneous_sidebar')

        if st.button("Gerar Sistema", key='generate_system_button'):
            config = {
                'num_vars': gen_num_vars,
                'var_notation': 'xyz' if gen_var_notation == 'x, y, z...' else 'indexed',
                'use_fractions': gen_use_fractions,
                'use_complex': gen_use_complex,
                'use_roots': gen_use_roots,
                'solution_type': gen_solution_type,
                'is_homogeneous': is_homogeneous,
            }
            generated_system, solution = system_generator.generate_system(config)
            st.session_state.user_input = generated_system
            st.session_state.generated_info = {"text": generated_system, "solution": solution, "variables": system_generator.get_variables(config)}
            # Clear the processed flag to force re-analysis on the main page
            if 'system_processed' in st.session_state:
                del st.session_state['system_processed']

    with st.expander("⚙️ Configurações Globais", expanded=True):
        st.radio("Modo de Notação", ["LaTeX", "Texto Plano (Unicode)"], key='notation_mode')

    st.markdown("---")
    st.markdown("##### 1. Sistema de Equações")
    
    st.text_area("Digite as equações:", key='user_input', height=150)
    
    if st.button("Analisar Sistema", type="primary"):
        matrix, variables, error = parse_system_input(st.session_state.user_input)
        if error: 
            st.error(f"Erro de Análise: {error}")
            st.session_state.system_processed = False
        else:
            st.success("Sistema processado!")
            st.info(f"Variáveis detectadas: {', '.join(variables)}")
            st.session_state.update(system_processed=True, current_matrix=matrix, variables=variables, matrix_history=[matrix], history=[], redo_history=[], redo_matrix_history=[], run_id=st.session_state.run_id + 1)
            # Clear generated info if a new system is analyzed
            if 'generated_info' in st.session_state:
                del st.session_state['generated_info']

    st.markdown("---")
    st.markdown("##### 2. Sistema Inicial")
    if st.session_state.get('system_processed'):
        if st.session_state.get('notation_mode', 'LaTeX') == 'LaTeX':
            st.latex(to_latex_system(st.session_state.matrix_history[0], st.session_state.variables))
        else:
            st.code(to_plain_text_system(st.session_state.matrix_history[0], st.session_state.variables))

    st.markdown("---")
    st.markdown("##### Controles")
    c1, c2, c3, c4, c5 = st.columns(5)
    with c1: st.button("⬅️", on_click=undo_last_operation, help="Desfazer a última operação", disabled=len(st.session_state.get('history', [])) == 0, use_container_width=True)
    with c2: st.button("➡️", on_click=redo_next_operation, help="Refazer a última operação desfeita", disabled=len(st.session_state.get('redo_history', [])) == 0, use_container_width=True)
    with c3: st.button("⏭️", on_click=apply_next_step, help="Aplica o próximo passo lógico para resolver o sistema.", use_container_width=True, disabled=not st.session_state.get('system_processed'))
    with c4: st.button("⏬", on_click=solve_to_ref, help="Resolve o sistema automaticamente para a Forma Escalonada (REF).", use_container_width=True, disabled=not st.session_state.get('system_processed'))
    with c5: st.button("⚡", on_click=solve_all_steps, help="Resolve o sistema completamente para a Forma Escalonada Reduzida (RREF).", use_container_width=True, disabled=not st.session_state.get('system_processed'))
    
    st.markdown("---")
    st.markdown("##### Histórico de Operações")
    if st.session_state.get('history'):
        for i, op in reversed(list(enumerate(st.session_state.history))):
            op_text = op[1] if st.session_state.get('notation_mode', 'LaTeX') == 'LaTeX' else op[0]
            if st.session_state.get('notation_mode', 'LaTeX') == 'LaTeX':
                st.latex(rf"{i + 1}: {op_text}")
            else:
                st.markdown(f"`{i + 1}: {op_text}`")
    else: st.info("Nenhuma operação realizada.")
    
    if st.session_state.get('generated_info') and st.session_state.user_input == st.session_state.generated_info['text']:
        st.markdown("---")
        st.markdown("##### Solução do Sistema Gerado")
        solution = st.session_state.generated_info.get('solution')
        variables = st.session_state.generated_info.get('variables', [])
        if solution is not None:
            for i, var in enumerate(variables):
                st.latex(f"{to_subscript(var)} = {latex(solution[i,0])}")
        else:
            st.info("O sistema gerado é impossível.")

st.header("Workbench de Sistemas Lineares")
if st.session_state.get('system_processed'):
    matrix = st.session_state.current_matrix
    variables = st.session_state.variables
    num_rows, num_vars = matrix.shape[0], len(variables)
    show_latex = st.session_state.notation_mode == 'LaTeX'
    pivot_info = get_pivot_info(matrix)

    st.markdown("##### Matriz Interativa")
    col_widths = [0.5, 2.5, 3.5] + [1.0] * num_vars + [1.0, 0.8]
    header_titles = ["L", "Operação", "Parâmetros"] + [to_subscript(v) for v in variables] + ["b", "Análise"]
    
    header_cols = st.columns(col_widths)
    for col, title in zip(header_cols, header_titles):
        col.markdown(f"<p style='text-align: center;'><strong>{title}</strong></p>", unsafe_allow_html=True)
    st.markdown("<hr style='margin:0; padding:0; margin-bottom: 5px;'/>", unsafe_allow_html=True)

    for i in range(num_rows):
        row_cols = st.columns(col_widths)
        row_cols[0].markdown(f"**$L_{{{i+1}}}$**")
        
        op_type = row_cols[1].selectbox("Op", ["[Nenhuma]", "Troca de Linhas", "Multiplicação por Escalar", "Combinação Linear"], key=f"op_{i}_{st.session_state.run_id}", label_visibility="collapsed")
        with row_cols[2]:
            other_rows = [r for r in range(1, num_rows + 1) if r != i + 1]
            if op_type == "Troca de Linhas":
                p_cols = st.columns([4,1]);
                r2 = p_cols[0].selectbox("com", other_rows, key=f"swap_{i}_{st.session_state.run_id}", format_func=lambda x: f"L{number_to_subscript(x)}", label_visibility="collapsed");
                p_cols[1].button("✅", key=f"apply_swap_{i}_{st.session_state.run_id}", on_click=apply_operation, args=(op_type,), kwargs={"r1":i+1, "r2":r2})
            elif op_type == "Multiplicação por Escalar":
                p_cols = st.columns([4,1]);
                alpha = p_cols[0].text_input("α", "1", key=f"mult_{i}_{st.session_state.run_id}", label_visibility="collapsed");
                p_cols[1].button("✅", key=f"apply_mult_{i}_{st.session_state.run_id}", on_click=apply_operation, args=(op_type,), kwargs={"r":i+1, "alpha":alpha})
            elif op_type == "Combinação Linear":
                p_cols = st.columns([2,2,1]);
                lambda_val = p_cols[0].text_input("λ", "1", key=f"sum_{i}_{st.session_state.run_id}", label_visibility="collapsed");
                rj = p_cols[1].selectbox("+λ⋅", other_rows, key=f"sum_j_{i}_{st.session_state.run_id}", format_func=lambda x: f"L{number_to_subscript(x)}", label_visibility="collapsed");
                p_cols[2].button("✅", key=f"apply_sum_{i}_{st.session_state.run_id}", on_click=apply_operation, args=(op_type,), kwargs={"ri":i+1, "rj":rj, "factor_lambda":lambda_val})

        current_row = list(matrix.row(i))
        for j, elem in enumerate(current_row):
            col = row_cols[j + 3]
            if show_latex: col.latex(latex(elem))
            else: col.markdown(f"<p style='text-align: center; font-family: monospace; height: 2rem; display:flex; align-items:center; justify-content:center;'>{format_plain_text(elem)}</p>", unsafe_allow_html=True)
        
        analysis_col = row_cols[num_vars + 4]
        with analysis_col.popover("🔍", use_container_width=False):
            st.markdown(f"**Análise da Linha $L_{{{i+1}}}$**")
            row_analysis = analyze_row_factors(matrix, i, pivot_info)
            is_absurd = row_analysis['tipo'][0] == '❗️'
            st.markdown("**Equação:**")
            if show_latex: st.latex(to_latex_equation(current_row, variables, is_absurd))
            else: st.code(to_plain_text_equation(current_row, variables, is_absurd))
            st.markdown("**Status:**")
            c1, c2, c3 = st.columns(3)
            with c1:
                st.markdown("Tipo")
                render_badge(row_analysis['tipo'][0], row_analysis['tipo'][1])
            with c2:
                st.markdown("Unitário?")
                render_badge(row_analysis['unitario'][0], row_analysis['unitario'][1])
            with c3:
                st.markdown("Reduzida?")
                render_badge(row_analysis['reduzida'][0], row_analysis['reduzida'][1])
            st.markdown("---")
            st.markdown("- ✅: Bom, ⚠️: Atenção, ❌: Ruim, ❗️: Contradição, ✳️: Linha Nula")

    st.markdown("<hr style='margin:0; padding:0; margin-top: 5px;'/>", unsafe_allow_html=True)
    
    st.markdown("##### Progresso para Forma Escalonada Reduzida (RREF)")
    progress = calculate_rref_progress(matrix, variables)
    st.progress(progress)

    with st.expander("Status do Sistema e Legenda", expanded=True):
        st.markdown("**Análise da Linha (no pop-over '🔍')**")
        st.markdown("- **Tipo**: ✅ Linha de Pivô | ❗️ Contradição (0=k) | ✳️ Nula (0=0)")
        st.markdown("- **Unit?**: ✅ Pivô é 1 | ⚠️ Pivô não é 1")
        st.markdown("- **Reduz?**: ✅ Linha Reduzida (só o pivô não é nulo) | ❌ Linha Não-Reduzida")
        
        st.markdown("---")
        st.markdown("**Status Global da Matriz**")
        if is_reduced_row_echelon(matrix): 
            st.success("🎉 Matriz na Forma Escalonada Reduzida (RREF)")
        elif is_row_echelon(matrix):
            st.info("ℹ️ Matriz na Forma Escalonada (REF)")
        else:
            st.warning("Matriz ainda não está escalonada.")
        
        classification, justification = classify_system(matrix, len(variables))
        st.info(f"**{classification}** ({justification})")

    with st.expander("Ver Resolução Completa (para copiar)"):
        log_items = generate_solution_log_items()
        for item_type, content in log_items:
            if item_type and content:
                if item_type == 'header': st.markdown(content)
                elif item_type == 'matrix_str': st.code(content, language=None)
else:
    st.info("Digite um sistema de equações na barra lateral para começar.")