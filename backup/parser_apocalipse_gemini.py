import sympy
import re

def parse_to_sympy(s):
    """Converte uma string para uma expressão SymPy, lidando com 'i' para complexos."""
    return sympy.sympify(s, locals={'i': sympy.I})

def find_symbols_in_equations(equations_str):
    """Encontra todos os símbolos (potenciais variáveis) nas equações."""
    try:
        # Encontra todas as sequências de letras seguidas opcionalmente por números
        # Ex: x, y, z, x1, var2, etc. Ignora 'i' se for um número complexo.
        raw_symbols = set(re.findall(r'\b([a-zA-Z][a-zA-Z0-9]*)\b', equations_str))
        
        # Filtra palavras-chave ou funções conhecidas se necessário
        # (deixado em branco por enquanto, mas poderia incluir 'sin', 'cos', 'exp')
        keywords = {'i'} # 'i' é tratado como a unidade imaginária
        
        symbols = {s for s in raw_symbols if s not in keywords}
        
        if not symbols:
            return set(), "Nenhum símbolo de variável foi encontrado. As variáveis devem começar com uma letra."
            
        return symbols, None
    except Exception as e:
        return None, f"Erro ao analisar símbolos: {e}"


def create_system_matrix(equations_str, variable_names):
    """
    Cria a matriz aumentada e uma lista ordenada de variáveis a partir de uma string de equações.
    """
    if not equations_str.strip():
        return None, None, "A entrada de equações está vazia."
    if not variable_names:
        return None, None, "A lista de nomes de variáveis não pode ser vazia."

    try:
        # Ordena as variáveis para garantir uma ordem consistente das colunas
        # Ex: ['x', 'y', 'z'] ou ['x1', 'x2', 'x3']
        variables = sorted([sympy.Symbol(v) for v in variable_names], key=str)
        
        lines = equations_str.strip().split('\n')
        num_equations = len(lines)
        num_vars = len(variables)
        
        A = sympy.zeros(num_equations, num_vars)
        b = sympy.zeros(num_equations, 1)

        for i, line in enumerate(lines):
            line = line.strip()
            if not line:
                # Se a linha estiver vazia, preenche com zeros
                for j in range(num_vars):
                    A[i, j] = 0
                b[i, 0] = 0
                continue

            if '=' not in line:
                return None, None, f"A equação na linha {i+1} ('{line}') não contém um sinal de '='."

            lhs_str, rhs_str = line.split('=', 1)
            
            # Substitui notação implícita tipo '2x' por '2*x' para o parser do SymPy
            # Isso é feito de forma cuidadosa para não afetar nomes de variáveis como 'x1'
            for var in variables:
                var_str = str(var)
                # Regex para encontrar um número (inteiro ou decimal) seguido diretamente pelo nome da variável
                # Lookbehind '(?<=\d)' e lookahead '(?=[a-zA-Z])' garantem que não estamos no meio de um número ou palavra
                lhs_str = re.sub(f'(\d+)(\*{{{var_str}}})', f'\\1*{{{var_str}}}', lhs_str) # Garante que 2*x1 não se torne 2**x1
                lhs_str = re.sub(f'(\d+)({var_str})', f'\\1*\\2', lhs_str)


            lhs_expr = parse_to_sympy(lhs_str.strip())
            rhs_expr = parse_to_sympy(rhs_str.strip())

            # Expande a expressão para separar os termos (ex: 2*x + 3*y)
            lhs_expanded = sympy.expand(lhs_expr)
            
            # Isola o termo constante no lado esquerdo
            constant_term = lhs_expanded.as_coeff_add()[0]
            if not constant_term.has(*variables):
                rhs_expr -= constant_term
                lhs_expanded -= constant_term

            # Extrai os coeficientes de cada variável
            for j, var in enumerate(variables):
                A[i, j] = lhs_expanded.coeff(var)

            b[i, 0] = rhs_expr

        # Constrói a matriz aumentada
        augmented_matrix = A.row_join(b)
        
        # Retorna os nomes das variáveis como strings
        variable_names_str = [str(v) for v in variables]
        
        return augmented_matrix, variable_names_str, None

    except Exception as e:
        return None, None, f"Erro ao criar a matriz do sistema: {e}"
