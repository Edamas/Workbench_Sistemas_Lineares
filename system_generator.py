"""
This module provides functions to generate random systems of linear equations
based on a set of user-defined constraints.
"""
import random
import sympy
from sympy import sqrt, I, Rational
import re

def _get_normal_randint(min_val, max_val, mu=0, sigma=1.5):
    """
    Generates a random integer following a normal distribution, clamped within a range.
    Values closer to mu are more probable.
    """
    while True:
        num = round(random.gauss(mu, sigma))
        if min_val <= num <= max_val:
            return num

def _get_random_value(config, is_for_solution=False):
    """
    Generates a single random value using a compositional probability model.
    """
    use_fractions = config.get('use_fractions', False)
    use_complex = config.get('use_complex', False)
    use_roots = config.get('use_roots', False)

    # 1. Base number: integer or fraction
    frac_chance = 0.3 if is_for_solution else 0.2 
    if use_fractions and random.random() < frac_chance:
        den = _get_normal_randint(-5, 5)
        while den == 0 or abs(den) == 1:
            den = _get_normal_randint(-5, 5)
        num = Rational(_get_normal_randint(-5, 5), den)
    else:
        num = sympy.Integer(_get_normal_randint(-5, 5))

    # 2. Composition with a root
    if use_roots and random.random() < 0.2: # 1/5 chance
        base = random.choice([2, 3, 5]) # Using simpler roots
        if num != 0:
            num *= sqrt(base)

    # 3. Composition with imaginary unit
    if use_complex and random.random() < 0.2: # 1/5 chance
        if num != 0:
            num *= I
            
    return num

def _get_random_coeff(config):
    """
    Generates a single random coefficient for the A matrix.
    It ensures the coefficient is not zero most of the time.
    """
    val = _get_random_value(config, is_for_solution=False)
    if val == 0 and random.random() > 0.25:
        return _get_random_coeff(config)
    return val

def _get_random_solution_value(config):
    """
    Generates a single random value for the solution vector 'x'.
    """
    return _get_random_value(config, is_for_solution=True)

def _custom_str(expr):
    """Custom string conversion to produce robust single-line output
    with prettier symbols for roots and complex numbers."""
    s = str(expr)
    s = re.sub(r'sqrt\((\d+)\)', r'√\1', s)
    s = s.replace('I', 'ⅈ')
    return s

def _format_term(coeff, var, is_first_term, config):
    """Formats a single term (coefficient * variable) into a string based on config."""
    if coeff == 0:
        return ""

    mul_symbol_option = config.get('mul_symbol', 'Nenhum (implícito)')
    mul_parentheses_enabled = config.get('mul_parentheses', False)

    mul_str = ""
    if mul_symbol_option == "Asterisco (*)":
        mul_str = "*"
    elif mul_symbol_option == "Ponto (⋅)":
        mul_str = "⋅"
    
    sign = ""
    if not is_first_term:
        # Simplified sign logic
        sign = " - " if (hasattr(coeff, 'is_real') and coeff.is_real and coeff < 0) else " + "
    elif (hasattr(coeff, 'is_real') and coeff.is_real and coeff < 0):
        sign = "-"

    coeff_abs = abs(coeff)
    
    if coeff_abs == 1:
        return f"{sign}{var}"

    coeff_str = _custom_str(coeff_abs)
    
    wrap_in_parens = False
    is_ambiguous_type = not coeff_abs.is_integer

    if mul_parentheses_enabled:
        wrap_in_parens = True
    elif is_ambiguous_type and mul_str == "":
        wrap_in_parens = True

    if wrap_in_parens:
        coeff_str = f"({coeff_str})"

    return f"{sign}{coeff_str}{mul_str}{var}"


def _format_system_to_string(matrix, variables, config):
    """Converts a SymPy Matrix into a string representation of the system, applying formatting config."""
    lines = []
    for i in range(matrix.rows):
        line_parts = []
        is_first_term = True
        for j, var in enumerate(variables):
            term = _format_term(matrix[i, j], var, is_first_term, config)
            if term:
                if is_first_term and term.startswith(" + "):
                    term = term[3:]
                line_parts.append(term)
                is_first_term = False
        
        lhs = "".join(line_parts).strip()
        if not lhs: 
            lhs = "0"
        elif lhs.startswith("+ "):
            lhs = lhs[2:]

        rhs = _custom_str(matrix[i, -1])
        lines.append(f"{lhs} = {rhs}")
    return "\n".join(lines)


def _generate_spd_system(config, variables):
    """Generates a 'Sistema Possível e Determinado' (SPD) and its solution."""
    num_vars = len(variables)
    A = sympy.zeros(num_vars, num_vars)
    while True:
        for i in range(num_vars):
            for j in range(num_vars):
                A[i, j] = _get_random_coeff(config)
        if A.det() != 0:
            break
            
    x = sympy.zeros(num_vars, 1)
    for i in range(num_vars):
        x[i, 0] = _get_random_solution_value(config)
    
    b = A * x

    if config.get('is_homogeneous', False):
        b = sympy.zeros(b.rows, 1)

    return A.row_join(b), x

def _generate_spi_system(config, variables):
    """Generates a 'Sistema Possível e Indeterminado' (SPI) and a particular solution."""
    num_vars = len(variables)
    if num_vars <= 1:
        return _generate_spd_system(config, variables)
    
    num_eqs = num_vars - random.randint(1, num_vars - 1)
    A = sympy.zeros(num_eqs, num_vars)
    for i in range(num_eqs):
        for j in range(num_vars):
            A[i, j] = _get_random_coeff(config)
            
    x = sympy.zeros(num_vars, 1)
    for i in range(num_vars):
        x[i, 0] = _get_random_solution_value(config)
    
    b = A * x

    if config.get('is_homogeneous', False):
        b = sympy.zeros(b.rows, 1)

    return A.row_join(b), x

def _generate_si_system(config, variables):
    """Generates a 'Sistema Impossível' (SI)."""
    num_vars = len(variables)
    num_eqs = num_vars
    
    A_consistent = sympy.zeros(num_eqs, num_vars)
    for i in range(num_eqs):
        for j in range(num_vars):
            A_consistent[i, j] = _get_random_coeff(config)

    x = sympy.zeros(num_vars, 1)
    for i in range(num_vars):
        x[i, 0] = _get_random_solution_value(config)
    b_consistent = A_consistent * x

    if config.get('is_homogeneous', False):
        return A_consistent.row_join(sympy.zeros(b_consistent.rows, 1)), x

    last_row_coeffs = sympy.zeros(1, num_vars)
    last_row_b = sympy.Integer(0)
    
    for _ in range(random.randint(1, 3)):
        row_idx = random.randint(0, num_eqs - 1)
        multiplier = _get_random_coeff(config)
        if multiplier == 0: continue
        
        last_row_coeffs += A_consistent.row(row_idx) * multiplier
        last_row_b += b_consistent.row(row_idx)[0] * multiplier

    conflict = _get_random_coeff(config)
    while conflict == 0:
        conflict = _get_random_coeff(config)
    last_row_b += conflict

    A = A_consistent.row_insert(num_eqs, last_row_coeffs)
    b = b_consistent.row_insert(num_eqs, sympy.Matrix([last_row_b]))
    
    return A.row_join(b), None

def get_variables(config):
    """Determines the list of variable names based on the configuration."""
    num_vars = config.get('num_vars', 3)
    var_notation = config.get('var_notation', 'indexed')
    if var_notation == 'xyz':
        return ['x', 'y', 'z', 'w', 't', 'p', 'q', 'r', 's', 'u'][:num_vars]
    else:
        return [f"x{i+1}" for i in range(num_vars)]

def generate_system(config):
    """
    Generates a system of linear equations and its solution based on the provided configuration.
    """
    solution = None
    variables = get_variables(config)
    solution_type = config.get('solution_type', 'Determinado')

    if solution_type == 'Determinado':
        system_matrix, solution = _generate_spd_system(config, variables)
    elif solution_type == 'Indeterminado':
        system_matrix, solution = _generate_spi_system(config, variables)
    elif solution_type == 'Impossível':
        system_matrix, solution = _generate_si_system(config, variables)
    else:
         return "# Tipo de solução não suportado", None

    return _format_system_to_string(system_matrix, variables, config), solution
