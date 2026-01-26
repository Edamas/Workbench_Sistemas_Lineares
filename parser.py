import sympy
from sympy import sympify, sqrt, I, Rational
import re
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, auto_number
from collections import OrderedDict
import numpy as np

def _pre_process_for_parse(text: str) -> str:
    """
    Prepara uma string de equação para o parser do SymPy, tornando a multiplicação implícita explícita.
    A ordem das substituições é importante.
    """
    # 1. Substitui símbolos unicode e normaliza
    text = text.replace('−', '-')
    text = text.replace('ⅈ', 'i') # Normaliza o 'i' unicode para 'i' padrão
    text = re.sub(r'√(\d+\.?\d*)', r'sqrt(\1)', text)
    
    # 2. Adiciona '*' para multiplicação implícita com a unidade imaginária 'i' ou 'I'
    # Ex: ix -> I*x, -ix1 -> -I*x1. A letra 'I' é usada para consistência.
    # A word boundary \b é crucial para não substituir 'six' por 's*I*x'
    text = re.sub(r'(-?)\b([iI])([a-zA-Z_][a-zA-Z0-9_]*)', r'\1I*\3', text)
    
    # 3. Adiciona '*' entre números e variáveis/parênteses
    # Ex: '2x' -> '2*x', '3(' -> '3*('
    text = re.sub(r'([\d\.]+)\s*([a-zA-Z\(])', r'\1*\2', text)
    
    # 4. Adiciona '*' após parênteses de fechamento
    # Ex: '(x+1)x' -> '(x+1)*x', '(x+y)(a-b)' -> '(x+y)*(a-b)'
    text = re.sub(r'(\))\s*([\da-zA-Z\(])', r'\1*\2', text)

    # 5. Adiciona '*' para casos como '2i' ou 'xi'
    # Ex: 'x i' -> 'x*I', 3i -> 3*I
    text = re.sub(r'([a-zA-Z\d\)])\s*\b([iI])\b', r'\1*I', text)

    return text

def parse_to_sympy(s: str):
    """
    Analisa uma string para uma expressão SymPy, assumindo que a multiplicação
    foi tornada explícita pelo pré-processador.
    """
    s_processed = _pre_process_for_parse(s)
    
    local_dict = {'I': sympy.I, 'i': sympy.I, 'sqrt': sympy.sqrt, 'raiz': sympy.sqrt} 
    
    # A multiplicação implícita é removida das transformações para evitar conflitos
    transformations = standard_transformations + (auto_number,)
    
    try:
        expr = parse_expr(s_processed, local_dict=local_dict, transformations=transformations, global_dict=None, evaluate=True)
        return expr
    except Exception as e:
        raise ValueError(f"Expressão '{s}' (processada como '{s_processed}') inválida. Erro: {e}")

def parse_system_input(text_input):
    """
    Analisa o texto de entrada do usuário em um sistema de equações lineares,
    identificando variáveis e extraindo coeficientes.
    """
    equations_str = [line.strip() for line in text_input.split('\n') if line.strip()]
    if not equations_str: return None, None, "Nenhuma equação encontrada."

    # --- FASE 1: Identificar todas as variáveis (símbolos) ---
    all_symbols_in_order = OrderedDict()
    detect_transformations = standard_transformations + (auto_number,)

    # Constrói uma única string de expressão para encontrar todas as variáveis de uma vez
    expression_parts = []
    for eq_str in equations_str:
        if '=' not in eq_str:
            return None, None, f"Erro na linha: '{eq_str}'. Falta o sinal '='."
        lhs, rhs = eq_str.split('=', 1)
        expression_parts.append(f"({lhs}) - ({rhs})")
    
    full_expression_string = " + ".join(expression_parts)
    full_processed_text = _pre_process_for_parse(full_expression_string)

    try:
        expr = parse_expr(full_processed_text, local_dict={'I': sympy.I, 'i': sympy.I, 'sqrt': sympy.sqrt}, transformations=detect_transformations, evaluate=False)
        
        # Use a natural sort key to handle variables like x1, x2, x10 correctly
        def natural_sort_key(s):
            return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', str(s))]

        sorted_symbols = sorted(list(expr.free_symbols), key=natural_sort_key)

        for s in sorted_symbols:
            if s != sympy.I and isinstance(s, sympy.Symbol): 
                if str(s) not in all_symbols_in_order:
                    all_symbols_in_order[str(s)] = s 
    except Exception as e:
        return None, None, f"Erro ao analisar o sistema para encontrar variáveis: {e}"
    
    variables = list(all_symbols_in_order.keys()) 
    if not variables:
        return None, None, "Nenhuma variável encontrada nas equações."
    
    sympy_variables = [sympy.Symbol(v) for v in variables]

    # --- FASE 2: Extrair Coeficientes e Termos Constantes ---
    matrix_np = np.empty((len(equations_str), len(variables) + 1), dtype=object) 

    for i, eq_str in enumerate(equations_str):
        lhs_str, rhs_str = eq_str.split('=', 1)
        
        try:
            lhs_expr_parsed = parse_to_sympy(lhs_str)
            rhs_expr_parsed = parse_to_sympy(rhs_str)
            
            # This is where the TypeError for tuples occurs
            equation_expr = lhs_expr_parsed - rhs_expr_parsed
            
            for j, var_symbol in enumerate(sympy_variables):
                matrix_np[i, j] = equation_expr.coeff(var_symbol)
            
            constant_term = equation_expr.as_independent(*sympy_variables)[0]
            matrix_np[i, -1] = -constant_term 

        except TypeError as e:
            if "unsupported operand type(s) for -" in str(e):
                error_msg = f"Na linha {i+1}: A notação de vetores/tuplas como `(a,b)` não é suportada para definir equações. Por favor, insira equações lineares padrão (ex: 3x + 2y = 5)."
                return None, None, error_msg
            else:
                # Re-raise other TypeErrors
                return None, None, f"Erro de tipo inesperado na linha {i+1}: {e}"
        except Exception as e:
            return None, None, f"Erro ao extrair coeficientes da equação '{eq_str}': {e}"

    # Garante que todas as entradas são objetos SymPy
    for r, c in np.ndindex(matrix_np.shape):
        if matrix_np[r, c] is None:
            matrix_np[r, c] = sympy.S.Zero
        else:
            matrix_np[r, c] = sympify(matrix_np[r, c]) 

    return sympy.Matrix(matrix_np), variables, None