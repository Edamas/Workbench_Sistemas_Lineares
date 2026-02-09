import sympy
from sympy import latex, Matrix

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

def generate_solution_log_items(matrix_history, op_history, user_input, variables):
    # This function is adapted from app.py and might need adjustments
    # to fit the data structures in app2.py.
    items = []
    if user_input:
        items.append(('header', "**Sistema de Equações Original:**"))
        items.append(('matrix_str', user_input))

    if matrix_history:
        items.append(('header', "**Matriz Aumentada Inicial [A|b]:**"))       
        items.append(('matrix_str', compact_format(matrix_history[0])))

    if not op_history:
        initial_matrix = matrix_history[0]
        if is_reduced_row_echelon(initial_matrix):
            items.append(('header', "**Observação:** Nenhuma operação foi necessária."))
    else:
        for i, op in enumerate(op_history):
            if op:
                items.append(('header', f"**Passo {i+1}:** `{op}`"))    
                items.append(('matrix_str', compact_format(matrix_history[i+1])))

    return items

def compact_format(matrix):
    formatted_rows = []
    for row_idx in range(matrix.rows):
        formatted_elements = [format_plain_text(elem) for elem in matrix.row(row_idx)]
        formatted_rows.append(f"[{', '.join(formatted_elements)}]")       
    return "\n".join(formatted_rows)

def format_plain_text(expr):
    if expr is None: return ""
    try:
        if expr.is_Integer: return str(expr)
        if expr.is_Float:
            val = float(expr)
            if val == round(val): return str(int(round(val)))       
            return f"{val:.3f}".rstrip('0').rstrip('.')
        if expr.is_Rational: return str(expr)
    except (AttributeError, TypeError): pass
    return sympy.pretty(expr, use_unicode=True)

def perform_subspace_axiom_test(axiom_name, dados_exercicio, counterexample_data=None):
    # This is a placeholder function. The actual implementation is missing.
    return None, "Teste automático ainda não implementado."
