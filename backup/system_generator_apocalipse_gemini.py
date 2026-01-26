import sympy
import random

def get_variables(config):
    """Gera uma lista de nomes de variáveis com base na configuração."""
    num_vars = config['num_vars']
    if config['var_notation'] == 'xyz':
        base_vars = ['x', 'y', 'z', 'w', 't', 'p', 'q', 'r', 's', 'u']
        return base_vars[:num_vars]
    else: # 'indexed'
        return [f"x{i+1}" for i in range(num_vars)]

def format_term(coeff, var_name, config):
    """Formata um único termo (ex: '2x', '-y', '3.5*z') com base na configuração."""
    if coeff == 0:
        return ""
    
    # Formatação do coeficiente
    coeff_str = ""
    if coeff == 1:
        coeff_str = ""
    elif coeff == -1:
        coeff_str = "-"
    else:
        # TODO: Implementar formatação para frações, raízes, etc.
        coeff_str = str(coeff)

    # Símbolo de multiplicação
    mul_symbol = ""
    if config['mul_symbol'] == "Asterisco (*)":
        mul_symbol = "*"
    elif config['mul_symbol'] == "Ponto (⋅)":
        mul_symbol = "⋅"
    
    # Parênteses no multiplicador
    if config['mul_parentheses'] and coeff not in [1, -1]:
        coeff_str = f"({coeff_str})"

    # Junta tudo
    if not coeff_str and mul_symbol: # Caso de coeficiente 1 com símbolo de mult
        return f"{var_name}"
    if coeff_str == "-":
        return f"-{var_name}"
    
    # Se não há coeficiente (é 1), não precisa de símbolo de multiplicação
    if not coeff_str:
        return var_name

    return f"{coeff_str}{mul_symbol}{var_name}"

def generate_equation(variables, solution_vector, config, is_homogeneous=False):
    """Gera uma única equação com base nos coeficientes e configuração."""
    terms = []
    # Gera coeficientes aleatórios para as variáveis
    coeffs = [random.randint(-10, 10) for _ in variables]
    
    # Calcula o lado direito (b) para que a solução seja válida
    rhs = sum(c * s for c, s in zip(coeffs, solution_vector))
    if is_homogeneous:
        rhs = 0

    # Formata cada termo
    for coeff, var in zip(coeffs, variables):
        term = format_term(coeff, var, config)
        if term:
            terms.append(term)
    
    # Junta os termos, cuidando dos sinais
    lhs_str = ""
    if not terms:
        lhs_str = "0"
    else:
        lhs_str = terms[0]
        for term in terms[1:]:
            if term.startswith('-'):
                lhs_str += f" {term}"
            else:
                lhs_str += f" + {term}"
                
    return f"{lhs_str} = {rhs}"


def generate_system(config):
    """
    Gera um sistema de equações lineares com base na configuração fornecida.
    """
    num_vars = config['num_vars']
    variables = get_variables(config)
    solution_type = config['solution_type']
    
    # 1. Gerar uma solução aleatória (ponto no R^n)
    # TODO: Lidar com frações, complexos, raízes na geração da solução
    solution = [random.randint(-5, 5) for _ in range(num_vars)]
    solution_matrix = sympy.Matrix(solution)

    equations = []
    num_equations = num_vars # Começa com um sistema quadrado
    
    if solution_type == "Indeterminado":
        # Para ser indeterminado, o número de equações independentes < num_vars
        num_equations = random.randint(1, num_vars - 1) if num_vars > 1 else 1
    elif solution_type == "Impossível":
        # Para ser impossível, podemos criar uma contradição
        # Gera n-1 equações consistentes e 1 inconsistente
        num_equations = num_vars

    # Gera as equações consistentes
    for _ in range(num_equations):
        eq = generate_equation(variables, solution, config, is_homogeneous=config['is_homogeneous'])
        equations.append(eq)
        
    # Adiciona manipulações para SPI e SI
    if solution_type == "Indeterminado" and num_vars > 1:
        # Adiciona uma ou mais equações que são combinações lineares das outras
        num_dependent_eqs = random.randint(1, num_vars - num_equations)
        # TODO: Implementar a lógica para criar combinações lineares
        pass

    elif solution_type == "Impossível" and num_vars > 0:
        # Pega a última equação e a torna inconsistente
        if equations:
            last_eq = equations[-1]
            lhs, rhs = last_eq.split('=')
            try:
                # Adiciona um valor aleatório não-zero ao lado direito
                rhs_val = int(rhs.strip())
                offset = random.choice([i for i in range(-5, 6) if i != 0])
                new_rhs = rhs_val + offset
                equations[-1] = f"{lhs.strip()} = {new_rhs}"
            except ValueError:
                # Se o RHS não for um inteiro simples, apenas anexa algo
                equations[-1] = f"{lhs.strip()} = {rhs.strip()} + 1"
        else: # Caso não tenha equações ainda (num_vars=1, impossível)
             eq = generate_equation(variables, solution, config)
             lhs, rhs = eq.split('=')
             equations.append(f"{lhs.strip()} = {int(rhs.strip()) + 1}")


    system_str = "\n".join(equations)
    
    # Para SPD e SPI, a solução é a que foi gerada. Para SI, não há solução.
    final_solution = solution_matrix if solution_type != "Impossível" else None
    
    return system_str, final_solution
