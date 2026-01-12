import streamlit as st
import sympy
from sympy import latex, Matrix
import system_generator
import numpy as np # Keep numpy for is_close in format_plain_text
import re # Keep re for to_subscript

# Import parsing functions from parser.py
from parser import parse_to_sympy, parse_system_input

# --- Configuração da Página ---
st.set_page_config(layout="wide", page_title="Workbench de Sistemas Lineares")

# --- Lógica de Análise e Matemática (Backend com SymPy) ---

def get_pivot_indices(matrix):
    pivots = {}
    if not hasattr(matrix, 'shape'): return pivots
    num_rows, num_cols = matrix.shape
    for r_idx in range(num_rows):
        row_coeffs = matrix.row(r_idx)
        for c_idx in range(num_cols - 1):
            if row_coeffs[c_idx] != 0:
                pivots[r_idx] = c_idx
                break
    return pivots

def is_row_echelon(matrix):
    pivots = get_pivot_indices(matrix)
    last_pivot_col = -1
    for r_idx in sorted(pivots.keys()):
        pivot_col = pivots[r_idx]
        if pivot_col <= last_pivot_col: return False
        for below_r_idx in range(r_idx + 1, matrix.shape[0]):
            if matrix[below_r_idx, pivot_col] != 0: return False
        last_pivot_col = pivot_col
    return True

def is_reduced_row_echelon(matrix):
    if not is_row_echelon(matrix): return False
    pivots = get_pivot_indices(matrix)
    for r_idx, pivot_col in pivots.items():
        if matrix[r_idx, pivot_col] != 1: return False
        for i in range(matrix.shape[0]):
            if i != r_idx and matrix[i, pivot_col] != 0: return False
    return True

def analyze_matrix_columns(matrix, variables):
    pivots = get_pivot_indices(matrix)
    pivot_cols = {pivots[r] for r in pivots}
    col_statuses = []
    for j in range(len(variables)):
        if j in pivot_cols:
            col_statuses.append(("Pivô", "success"))
        else:
            col_statuses.append(("Livre", "info"))
    return col_statuses

def analyze_row_sympy(matrix, row_index, variables):
    statuses = {
        'pivot': ('—', 'normal'),
        'col_exclusive': ('—', 'normal'),
        'row_exclusive': ('—', 'normal')
    }
    pivots = get_pivot_indices(matrix)
    num_vars = len(variables)

    if row_index in pivots:
        pivot_col = pivots[row_index]
        pivot_value = matrix[row_index, pivot_col]

        # 1. Pivot Status
        if pivot_value == 1:
            statuses['pivot'] = ("Unitário", "success")
        else:
            statuses['pivot'] = ("Não Unitário", "warning")

        # 2. Exclusive Column Status
        is_col_cleared = all(matrix[r, pivot_col] == 0 for r in range(matrix.shape[0]) if r != row_index)
        statuses['col_exclusive'] = ("Sim", "success") if is_col_cleared else ("Não", "error")

        # 3. Exclusive Row Status
        is_row_reduced = all(matrix[row_index, c] == 0 for c in range(num_vars) if c != pivot_col)
        statuses['row_exclusive'] = ("Sim", "success") if is_row_reduced else ("Não", "info")
        
    else: # No pivot in this row
        if all(c == 0 for c in matrix[row_index, :-1]):
            if matrix[row_index, -1] != 0: 
                statuses['pivot'] = ("Contradição", "error")
            else: 
                statuses['pivot'] = ("Nula", "info")
    
    return statuses

def calculate_rref_progress(matrix, variables):
    num_rows = matrix.shape[0]
    if num_rows == 0: return 0.0
    
    rref_compliant_rows = 0
    for i in range(num_rows):
        analysis = analyze_row_sympy(matrix, i, variables)
        pivot_status = analysis['pivot'][0]
        col_status = analysis['col_exclusive'][0]
        row_status = analysis['row_exclusive'][0]
        
        if pivot_status == "Nula":
            rref_compliant_rows += 1
        elif pivot_status == "Unitário" and col_status == "Sim" and row_status == "Sim":
            rref_compliant_rows += 1
            
    return rref_compliant_rows / num_rows


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

def to_subscript(s):
    match = re.match(r'([a-zA-Z])(\d+)', s)
    if not match: return s
    name, digits = match.groups()
    subscript_digits = "".join(['₀₁₂₃₄₅₆₇₈₉'[int(d)] for d in digits])
    return name + subscript_digits

def number_to_subscript(n):
    """Converts a number to a string of subscript digits."""
    return "".join(['₀₁₂₃₄₅₆₇₈₉'[int(d)] for d in str(n)])

def format_plain_text(expr):
    """Formats a SymPy expression for plain text, handling rounding and unicode."""
    if expr is None: return ""
    try:
        if expr.is_Integer:
            return str(expr)
        if expr.is_Float:
            val = float(expr)
            if np.isclose(val, round(val)):
                return str(int(round(val)))
            return f"{val:.3f}".rstrip('0').rstrip('.')
        if expr.is_number: # Handles Rational, Complex, etc.
            # If complex number with float parts, format it slightly
            if expr.is_complex and (expr.as_real_imag()[0].is_Float or expr.as_real_imag()[1].is_Float):
                real_part = float(expr.as_real_imag()[0])
                imag_part = float(expr.as_real_imag()[1])
                # Format real and imag parts recursively for consistent rounding
                formatted_real = format_plain_text(real_part)
                formatted_imag = format_plain_text(abs(imag_part))
                if np.isclose(imag_part, 0): return formatted_real
                if np.isclose(real_part, 0): return f"{formatted_imag}i" if imag_part > 0 else f"-{{formatted_imag}}i"
                return f"({{formatted_real}}{'+' if imag_part > 0 else '-'}{formatted_imag}i)"
            return str(expr) # Default string representation for other numeric types
    except (AttributeError, TypeError):
        pass
    
    pretty_str = sympy.pretty(expr, use_unicode=True)
    if "\n" in pretty_str:
        return str(expr)
    return pretty_str

def to_latex_equation(row_vector, variables):
    line = []
    for j, var_name in enumerate(variables):
        coeff = row_vector[j]
        if coeff != 0:
            term_latex = latex(coeff * sympy.Symbol(var_name), mul_symbol=' ')
            if coeff == 1: term_latex = latex(sympy.Symbol(var_name))
            elif coeff == -1: term_latex = f"-{latex(sympy.Symbol(var_name))}"
            line.append(term_latex if not line else (f"+ {term_latex}" if not term_latex.startswith('-') else term_latex))
    lhs = " ".join(line).replace("+ -", "- ").strip() or "0"
    return f"{lhs} = {latex(row_vector[-1])}"

def to_latex_system(matrix, variables):
    # Add S = before the system and use correct raw strings for LaTeX
    begin_system = r"S = \begin{cases}"
    end_system = r"\end{cases}"
    content = r" \\ ".join([to_latex_equation(list(matrix.row(i)), variables) for i in range(matrix.shape[0])])
    return begin_system + content + end_system

def apply_operation(op_type, **kwargs):
    st.session_state.redo_history = []
    st.session_state.redo_matrix_history = []
    current_step = len(st.session_state.history)
    if current_step < len(st.session_state.matrix_history) -1:
        st.session_state.matrix_history = st.session_state.matrix_history[:current_step + 1]
        st.session_state.history = st.session_state.history[:current_step]
        st.sidebar.warning("Histórico futuro descartado.")
    
    matrix = st.session_state.current_matrix.copy()
    variables = st.session_state.variables
    try:
        plain_desc = ""
        latex_desc = ""

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
            for j in range(matrix.cols): matrix[ri, j] = sympy.simplify(matrix[ri, j] + lambda_val * matrix[rj, j])
            plain_desc = f"L{ri+1} <- L{ri+1} + ({lambda_str}) * L{rj+1}"
            latex_desc = rf"L_{{{ri+1}}} \leftarrow L_{{{ri+1}}} + ({latex(lambda_val)}) \cdot L_{{{rj+1}}}"
        else:
            return

        # Add Op Type abbreviation
        op_type_abbr = {
            "Troca de Linhas": "Troca",
            "Multiplicação por Escalar": "Mult.",
            "Combinação Linear": "Comb. Lin."
        }.get(op_type, op_type)
        plain_desc += f" [{op_type_abbr}]"

        st.session_state.current_matrix = matrix

        # Check for echelon form and append to description
        if is_reduced_row_echelon(matrix):
            plain_desc += " [Matriz na Forma Escalonada Reduzida]"
            st.sidebar.success("🎉 Matriz na forma escalonada reduzida!")
        elif is_row_echelon(matrix):
            plain_desc += " [Matriz na Forma Escalonada]"
            st.sidebar.info("ℹ️ Matriz na forma escalonada.")

        st.session_state.history.append((plain_desc, latex_desc))
        st.session_state.matrix_history.append(matrix)
        st.sidebar.success("Operação aplicada!")
        st.sidebar.latex(latex_desc)
    except Exception as e: st.sidebar.error(f"Erro na operação: {e}")

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

def replay_to_step(step_index):
    if 0 <= step_index < len(st.session_state.matrix_history):
        while len(st.session_state.matrix_history) > step_index + 1:
            st.session_state.redo_matrix_history.append(st.session_state.matrix_history.pop())
            st.session_state.redo_history.append(st.session_state.history.pop())
        st.session_state.current_matrix = st.session_state.matrix_history[step_index]
        st.sidebar.info(f"Visualizando estado após a operação {step_index}.")

def render_badge(badge_type, text, icon_only=False):
    icon_map = {"success": "✅", "warning": "⚠️", "error": "❌", "info": "ℹ️"}
    icon = icon_map.get(badge_type)
    
    if icon_only:
        st.markdown(f"<div style='text-align: center;'>{icon}</div>", unsafe_allow_html=True)
    elif badge_type == "success":
        st.success(text, icon=icon)
    elif badge_type == "warning":
        st.warning(text, icon=icon)
    elif badge_type == "error":
        st.error(text, icon=icon)
    elif badge_type == "info":
        st.info(text, icon=icon)
    else:
        st.markdown(text)

def to_plain_text_equation(row_vector, variables):
    line = []
    for j, var_name in enumerate(variables):
        coeff, var_symbol = row_vector[j], sympy.Symbol(var_name)
        if coeff != 0:
            if coeff == 1: 
                term_text = to_subscript(var_name)
            elif coeff == -1: 
                term_text = f"-{to_subscript(var_name)}"
            else: 
                term_text = f"{format_plain_text(coeff)}*{to_subscript(var_name)}"
            
            if not line:
                line.append(term_text)
            elif term_text.startswith("-"):
                 line.append(f" {term_text}")
            else:
                 line.append(f" + {term_text}")

    lhs = "".join(line).strip() or "0"
    return f"{lhs} = {format_plain_text(row_vector[-1])}"

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
    
    original_input_str = st.session_state.get('original_input', '')
    if original_input_str:
        items.append(('header', "**Sistema de Equações Original:**"))
        items.append(('matrix_str', original_input_str))

    if st.session_state.get('matrix_history'):
        items.append(('header', "**Matriz Aumentada Inicial [A|b]:**"))
        items.append(('matrix_str', compact_format(st.session_state.matrix_history[0])))

    for i, (plain_op, _) in enumerate(st.session_state.history):
        if plain_op: 
            items.append(('header', f"**Passo {i+1}:** `{plain_op}`"))
            items.append(('matrix_str', compact_format(st.session_state.matrix_history[i+1])))
    
    items.append(('header', "---\n\n**ANÁLISE FINAL DO SISTEMA**\n"))
    final_matrix = st.session_state.current_matrix
    variables = st.session_state.variables
    num_vars = len(variables)
    
    if is_reduced_row_echelon(final_matrix):
        final_matrix_label = "**1. Matriz Final na Forma Escalonada Reduzida:**"
    else:
        final_matrix_label = "**1. Matriz Final na Forma Escalonada:**"
    items.append(('header', final_matrix_label))
    items.append(('matrix_str', compact_format(final_matrix)))

    rank_A, rank_Ab = final_matrix[:,:-1].rank(), final_matrix.rank()
    analysis_text = f"**2. Análise dos Postos:**\n- Posto da matriz de coeficientes: **posto(A) = {rank_A}**\n- Posto da matriz aumentada: **posto(A|b) = {rank_Ab}**\n- Número de incógnitas: **n = {num_vars}**\n\n**3. Conclusão (Teorema de Rouché-Capelli):**"
    items.append(('header', analysis_text))
    if rank_A != rank_Ab:
        items.append(('header', f"Como posto(A) ≠ posto(A|b), o sistema é **Impossível (SI)**.\n- **Conjunto Solução (S):** S = ∅"))
    else:
        if rank_A == num_vars:
            items.append(('header', f"Como posto(A) = posto(A|b) = n = {num_vars}, o sistema é **Possível e Determinado (SPD)**."))
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

# --- Interface Streamlit (Frontend) ---
if 'system_processed' not in st.session_state:
    st.session_state.system_processed = False

# Ensure history lists are initialized
if 'history' not in st.session_state:
    st.session_state.history = []
    st.session_state.matrix_history = []
    st.session_state.redo_history = []
    st.session_state.redo_matrix_history = []
if 'run_id' not in st.session_state:
    st.session_state.run_id = 0

with st.sidebar:
    st.title("Workbench de Sistemas Lineares")
    st.markdown("---")

    with st.expander("0. Gerador de Sistema", expanded=False):
        gen_num_vars = st.slider("Número de variáveis", 1, 10, 3)
        gen_var_notation = st.radio("Notação das variáveis", ["x, y, z...", "x₁, x₂, x₃..."], index=1)
        
        gen_mul_symbol = st.radio("Símbolo de Multiplicação", ["Nenhum (implícito)", "Asterisco (*)", "Ponto (⋅)"], index=0, help="Escolha como o multiplicador será exibido (ex: 2x, 2*x, 2⋅x).")
        gen_mul_parentheses = st.radio("Multiplicador entre Parênteses", ["Não", "Sim"], index=0, help="Ex: (2)x ou 2x. Apenas para multiplicadores diferentes de 1 ou -1.")

        c1, c2, c3 = st.columns(3)
        gen_use_fractions = c1.checkbox("Frações", value=False)
        gen_use_complex = c2.checkbox("Complexos")
        gen_use_roots = c3.checkbox("Raízes")
        
        gen_solution_type = st.selectbox("Tipo de Solução", ["Determinado", "Indeterminado", "Impossível"])
        is_homogeneous = st.checkbox("Sistema homogêneo", value=False)

        if st.button("Gerar Sistema"):
            config = {
                'num_vars': gen_num_vars,
                'var_notation': 'xyz' if gen_var_notation == 'x, y, z...' else 'indexed',
                'mul_symbol': gen_mul_symbol,
                'mul_parentheses': gen_mul_parentheses == 'Sim',
                'use_fractions': gen_use_fractions,
                'use_complex': gen_use_complex,
                'use_roots': gen_use_roots,
                'solution_type': gen_solution_type,
                'is_homogeneous': is_homogeneous,
            }
            generated_system, solution = system_generator.generate_system(config)
            st.session_state.user_input = generated_system
            st.session_state.generated_info = {"text": generated_system, "solution": solution, "variables": system_generator.get_variables(config)}

    st.markdown("---")

    st.markdown("##### 1. Sistema de Equações")
    if 'user_input' not in st.session_state:
        st.session_state.user_input = "x + 3y - 2z = 4 - 4i\n-ix + 2y + z = 8\nx + y - z = 1"

    user_input = st.text_area("Digite as equações:", key='user_input', height=150, help="Use '*' para multiplicação (ex: 3*x2) ou notação implícita (3x2). Use '^' ou '**' para potências.")
    
    if st.button("Analisar Sistema", type="primary"):
        # Increment run_id to force widget re-creation
        run_id = st.session_state.get('run_id', 0) + 1
        current_input = st.session_state.user_input
        
        # Preserve generated info if input hasn't changed
        gen_info = st.session_state.get('generated_info')
        if gen_info and gen_info.get('text') != current_input:
            gen_info = None

        st.session_state.clear()
        
        # Restore essential state
        st.session_state.user_input = current_input
        st.session_state.run_id = run_id
        if gen_info:
            st.session_state.generated_info = gen_info
        
        matrix, variables, error = parse_system_input(current_input)
        if error: 
            st.error(f"Erro de Análise: {error}")
            st.session_state.system_processed = False
        else:
            st.success("Sistema processado!")
            st.session_state.update(system_processed=True, current_matrix=matrix, variables=variables, matrix_history=[matrix], history=[], original_input=current_input, redo_history=[], redo_matrix_history=[])

    st.markdown("---")
    st.markdown("##### 2. Sistema Inicial")
    if st.session_state.get('system_processed'): st.latex(to_latex_system(st.session_state.matrix_history[0], st.session_state.variables))
    else: st.info("Aguardando análise do sistema.")
    st.markdown("---")
    st.markdown("##### 3. Histórico de Operações")
    undo_col, redo_col = st.columns(2)
    with undo_col:
        st.button("⬅️", on_click=undo_last_operation, help="Desfaz a última operação (Ctrl+Z)", disabled=len(st.session_state.get('history', [])) == 0)
    with redo_col:
        st.button("➡️", on_click=redo_next_operation, help="Refaz a última operação desfeita (Ctrl+Y)", disabled=len(st.session_state.get('redo_history', [])) == 0)

    if st.session_state.get('history'):
        for i, op in reversed(list(enumerate(st.session_state.history))):
            cols = st.columns([4,1]); cols[0].latex(rf"{i + 1}: {op[1]}")
            cols[1].button("Ver", key=f"replay_{i}", on_click=replay_to_step, args=(i + 1,), help=f"Voltar para o estado após a operação {i+1}")
    else: st.info("Nenhuma operação realizada.")
    st.markdown("---")

    if st.session_state.get('generated_info') and st.session_state.user_input == st.session_state.generated_info['text']:
        st.markdown("##### 4. Solução do Sistema Gerado")
        solution = st.session_state.generated_info.get('solution')
        variables = st.session_state.generated_info.get('variables', [])
        if solution is not None:
            for i, var in enumerate(variables):
                st.latex(f"{to_subscript(var)} = {latex(solution[i,0])}")
        else:
            st.info("O sistema gerado é impossível (não possui solução).")

st.header("Workbench de Sistemas Lineares")
if st.session_state.get('system_processed'):
    matrix, variables = st.session_state.current_matrix, st.session_state.variables
    num_rows = matrix.shape[0]

    # --- Progress Bar ---
    st.markdown("##### Progresso para Forma Escalonada Reduzida (RREF)")
    progress = calculate_rref_progress(matrix, variables)
    st.progress(progress)
    
    # --- Column Analysis Header ---
    st.markdown("###### Análise das Colunas")
    # Adicionada uma coluna extra para LaTeX
    cef_cols = st.columns([1, 2, 3] + [1] * len(variables) + [1, 3, 3, 0.8, 1.2, 1.2])
    cef_cols[0].markdown("**Status**")
    col_statuses = analyze_matrix_columns(matrix, variables)
    for i, (text, badge_type) in enumerate(col_statuses):
        with cef_cols[i + 3]:
            render_badge(badge_type, text)
    
    st.markdown("---")

    # --- Main Workbench Header ---
    display_variables = [to_subscript(v) for v in variables]
    header_cols = st.columns([1, 2, 3] + [1] * len(variables) + [1, 3, 3, 0.8, 1.2, 1.2])
    # Adicionado "Equação (LaTeX)" ao cabeçalho
    for col, title in zip(header_cols, ["Linha", "Operação", "Parâmetros"] + display_variables + ["b", "Equação", "Equação (LaTeX)", "Pivô", "Col. Excl?", "Linha Excl?"]): 
        col.markdown(f"**{title}**")
    st.markdown("<hr style='margin-top: 0; margin-bottom: 0;'/>", unsafe_allow_html=True)

    # --- Main Workbench Rows ---
    for i in range(num_rows):
        row_cols = st.columns([1, 2, 3] + [1] * len(variables) + [1, 3, 3, 0.8, 1.2, 1.2])
        row_cols[0].markdown(f"### $L_{{{i+1}}}$")
        op_type = row_cols[1].selectbox("Operação", ["[Nenhuma]", "Troca de Linhas", "Multiplicação por Escalar", "Combinação Linear"], key=f"op_type_{i}_{st.session_state.run_id}", label_visibility="collapsed", index=0)
        with row_cols[2]:
            other_rows = [r for r in range(1, num_rows + 1) if r != i + 1]
            if op_type == "Troca de Linhas":
                p_cols = st.columns([4,1]); 
                r2 = p_cols[0].selectbox("com L", other_rows, key=f"swap_j_{i}_{st.session_state.run_id}", format_func=lambda x: f"L{number_to_subscript(x)}", label_visibility="collapsed"); 
                p_cols[1].button("✅", key=f"apply_swap_{i}_{st.session_state.run_id}", on_click=apply_operation, args=(op_type,), kwargs={"r1":i+1, "r2":r2}, help="Aplicar Troca de Linhas")
            elif op_type == "Multiplicação por Escalar":
                p_cols = st.columns([4,1]); 
                alpha = p_cols[0].text_input("α", "1", key=f"mult_alpha_{i}_{st.session_state.run_id}", label_visibility="collapsed"); 
                p_cols[1].button("✅", key=f"apply_mult_{i}_{st.session_state.run_id}", on_click=apply_operation, args=(op_type,), kwargs={"r":i+1, "alpha":alpha}, help="Aplicar Multiplicação por Escalar")
            elif op_type == "Combinação Linear":
                p_cols = st.columns([2,2,1]); 
                lambda_val = p_cols[0].text_input("λ", "-1", key=f"sum_lambda_{i}_{st.session_state.run_id}", label_visibility="collapsed"); 
                rj = p_cols[1].selectbox("+λ⋅", other_rows, key=f"sum_j_{i}_{st.session_state.run_id}", format_func=lambda x: f"L{number_to_subscript(x)}", label_visibility="collapsed"); 
                p_cols[2].button("✅", key=f"apply_sum_{i}_{st.session_state.run_id}", on_click=apply_operation, args=(op_type,), kwargs={"ri":i+1, "rj":rj, "factor_lambda":lambda_val}, help="Aplicar Combinação Linear")
        
        current_matrix_row = list(matrix.row(i))
        for col_idx, elem in enumerate(current_matrix_row):
            row_cols[col_idx + 3].latex(latex(elem))
        
        # Coluna de Texto Plano
        row_cols[len(variables) + 4].code(to_plain_text_equation(current_matrix_row, variables))
        # Nova Coluna de LaTeX
        row_cols[len(variables) + 5].latex(to_latex_equation(current_matrix_row, variables))

        observations = analyze_row_sympy(matrix, i, variables)
        with row_cols[len(variables) + 6]:
            render_badge(observations['pivot'][1], "", icon_only=True)
        with row_cols[len(variables) + 7]:
            render_badge(observations['col_exclusive'][1], "", icon_only=True)
        with row_cols[len(variables) + 8]:
            render_badge(observations['row_exclusive'][1], "", icon_only=True)
        
        st.markdown("<hr style='margin-top: 0; margin-bottom: 0;'/>", unsafe_allow_html=True)
    
    # --- Legend and Status Expander ---
    with st.expander("Status do Sistema e Legenda", expanded=False):
        leg_col1, leg_col2, leg_col3 = st.columns([2, 2, 3])
        
        with leg_col1:
            st.markdown("""
            **Análise das Colunas**
            - ✅ **Pivô:** Coluna contém um pivô.
            - ℹ️ **Livre:** Coluna de variável livre.
            """, unsafe_allow_html=True)
        with leg_col2:
            st.markdown("""
            **Análise das Linhas**
            - **Pivô**: ✅ Unitário / ⚠️ Não-Unitário
            - **Col. Excl?**: ✅ Sim / ❌ Não
            - **Linha Excl?**: ✅ Sim / ℹ️ Não
            """, unsafe_allow_html=True)
        with leg_col3:
            st.markdown("**Status Global do Sistema**")
            if is_reduced_row_echelon(matrix): 
                st.success("🎉 Matriz na Forma Escalonada Reduzida (RREF)")
            elif is_row_echelon(matrix):
                st.info("ℹ️ Matriz na Forma Escalonada (REF)")
            else:
                st.warning("Matriz ainda não está escalonada.")
            
            classification, justification = classify_system(matrix, len(variables))
            st.info(f"**{classification}** ({justification})")

    # Display solution only if system is solved
    if is_reduced_row_echelon(matrix):
        st.markdown("---")
        st.header("Solução do Sistema")
        classification, _ = classify_system(matrix, len(variables))
        if classification == "Sistema Possível e Determinado (SPD)":
            solution, error = solve_echelon_spd(matrix, variables)
            if not error:
                has_complex = any(val.is_complex for val in solution.values() if val is not None and hasattr(val, 'is_complex'))
                field = "\\mathbb{C}" if has_complex else "\\mathbb{R}"
                ordered_solution_tuple = f"({', '.join([latex(solution.get(v, '')) for v in variables])})"
                solution_set_str = f"S = \\{{{ordered_solution_tuple}\\}} \\subset {field}^{{{len(variables)}}}"
                st.latex(solution_set_str)

        elif classification == "Sistema Possível e Indeterminado (SPI)":
            st.markdown("**Solução Geral:**")
            st.code(to_plain_text_system(matrix, variables))

    if st.session_state.get('history'):
        st.markdown("---")
        with st.expander("Ver Resolução Completa (para copiar)"):
            log_items = generate_solution_log_items()
            for item_type, content in log_items:
                if item_type and content:
                    if item_type == 'header': st.markdown(content)
                    elif item_type == 'matrix_str': st.code(content, language=None)
else:
    st.info("Digite um sistema de equações na barra lateral para começar.")