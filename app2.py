import streamlit as st
import json
import os
import sympy
from sympy import latex

# Importa os módulos necessários
import parser
import solver_utils

# Define o caminho para os arquivos de dados e resoluções
EXERCICIOS_FILE = 'exercicios.json'
RESOLUCOES_DIR = 'resolucoes'


# --- CONFIGURAÇÃO DA PÁGINA ---
st.set_page_config(
    page_title="Assistente de Álgebra Linear",
    page_icon="🤖",
    layout="wide"
)

# --- FUNÇÕES DE CARREGAMENTO DE DADOS ---

@st.cache_data
def carregar_exercicios():
    """Carrega a lista de exercícios do arquivo JSON."""
    if not os.path.exists(EXERCICIOS_FILE):
        st.error(f"Arquivo de exercícios não encontrado: {EXERCICIOS_FILE}")
        return []
    with open(EXERCICIOS_FILE, 'r', encoding='utf-8') as f:
        return json.load(f)

def carregar_resolucao(exercicio_id):
    """Carrega a resolução de um exercício, se existir."""
    resolucao_path = os.path.join(RESOLUCOES_DIR, f"{exercicio_id}.json")
    if os.path.exists(resolucao_path):
        with open(resolucao_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    return None

def salvar_resolucao(exercicio_id, resolucao_data):
    """Salva os dados da resolução em um arquivo JSON."""
    if not os.path.exists(RESOLUCOES_DIR):
        os.makedirs(RESOLUCOES_DIR)
    
    resolucao_path = os.path.join(RESOLUCOES_DIR, f"{exercicio_id}.json")
    with open(resolucao_path, 'w', encoding='utf-8') as f:
        json.dump(resolucao_data, f, indent=4, ensure_ascii=False)

# --- Funções de operação (globais, operam em st.session_state.solver_state) ---
def apply_op(op_type, source="manual", **kwargs):
    solver_state = st.session_state.solver_state # Obter o estado atual do solver
    if source == "manual":
        solver_state['redo_history'] = []
        solver_state['redo_matrix_history'] = []
    
    matrix = solver_state['current_matrix'].copy()
    try:
        plain_desc, latex_desc = "", ""
        if op_type == "Troca de Linhas":
            r1, r2 = kwargs['r1'] - 1, kwargs['r2'] - 1
            matrix.row_swap(r1, r2)
            plain_desc = f"L{r1+1} <-> L{r2+1}"
            latex_desc = f"L_{{{r1+1}}} \leftrightarrow L_{{{r2+1}}}"
        elif op_type == "Multiplicação por Escalar":
            r, alpha_str = kwargs['r'] - 1, kwargs['alpha']
            alpha = parser.parse_to_sympy(alpha_str)
            if alpha == 0: raise ValueError("O escalar (α) não pode ser zero.")
            for j in range(matrix.cols): matrix[r, j] = sympy.simplify(matrix[r, j] * alpha)
            plain_desc = f"L{r+1} <- ({alpha_str}) * L{r+1}"
            latex_desc = f"L_{{{r+1}}} \leftarrow ({latex(alpha)}) \cdot L_{{{r+1}}}"
        elif op_type == "Combinação Linear":
            ri, rj, lambda_str = kwargs['ri'] - 1, kwargs['rj'] - 1, kwargs['factor_lambda']
            lambda_val = parser.parse_to_sympy(lambda_str)
            if lambda_val == 0: raise ValueError("O fator (λ) não pode ser zero.")
            matrix.row_add(ri, rj, lambda_val)
            plain_desc = f"L{ri+1} <- L{ri+1} + ({lambda_str}) * L{rj+1}"
            latex_desc = f"L_{{{ri+1}}} \leftarrow L_{{{ri+1}}} + ({latex(lambda_val)}) \cdot L_{{{rj+1}}}"
        else: return

        solver_state['current_matrix'] = matrix
        solver_state['op_history'].append(latex_desc) # Usando latex_desc no histórico
        solver_state['matrix_history'].append(matrix)
        solver_state['run_id'] += 1 # Incrementa para forçar re-renderização de widgets

    except Exception as e:
        st.error(f"Erro na operação: {e}")
    st.rerun() # Força a atualização da UI

def undo_op():
    solver_state = st.session_state.solver_state # Obter o estado atual do solver
    if len(solver_state['op_history']) > 0:
        solver_state['redo_history'].append(solver_state['op_history'].pop())
        solver_state['redo_matrix_history'].append(solver_state['matrix_history'].pop())
        solver_state['current_matrix'] = solver_state['matrix_history'][-1]
        solver_state['run_id'] += 1 # Para forçar a atualização dos widgets
    st.rerun()

def redo_op():
    solver_state = st.session_state.solver_state # Obter o estado atual do solver
    if len(solver_state['redo_history']) > 0:
        solver_state['op_history'].append(solver_state['redo_history'].pop())
        solver_state['matrix_history'].append(solver_state['redo_matrix_history'].pop())
        solver_state['current_matrix'] = solver_state['matrix_history'][-1]
        solver_state['run_id'] += 1 # Para forçar a atualização dos widgets
    st.rerun()

def apply_next_auto_step():
    solver_state = st.session_state.solver_state # Obter o estado atual do solver
    op = solver_utils.find_next_step(solver_state['current_matrix'])
    if op:
        op_type = op.pop('op_type')
        apply_op(op_type, source="auto", **op)
    else:
        st.success("A matriz já está na RREF!")

def solve_to_ref_auto():
    solver_state = st.session_state.solver_state # Obter o estado atual do solver
    MAX_STEPS = 50 
    for _ in range(MAX_STEPS):
        op = solver_utils.find_next_ref_step(solver_state['current_matrix'])
        if op:
            op_type = op.pop('op_type')
            apply_op(op_type, source="auto", **op)
        else:
            st.success("Sistema escalonado (REF)!")
            return
    st.error(f"A resolução parou após {MAX_STEPS} passos.")

def solve_all_steps_auto():
    solver_state = st.session_state.solver_state # Obter o estado atual do solver
    MAX_STEPS = 50 
    for _ in range(MAX_STEPS):
        op = solver_utils.find_next_step(solver_state['current_matrix'])
        if op:
            op_type = op.pop('op_type')
            apply_op(op_type, source="auto", **op)
        else:
            st.success("Sistema resolvido para RREF!")
            return
    st.error(f"A resolução automática parou após {MAX_STEPS} passos.")


# --- RENDERIZAÇÃO DINÂMICA (DISPATCHER) ---

def renderizar_verificacao_subespaco(problema_data, problem_id="custom"):
    """Renderiza a interface para exercícios de verificação de subespaço com teste automatizado."""
    exercicio_id = problem_id # Usar problem_id para salvar e carregar
    dados_exercicio = problema_data # Usar problema_data diretamente
    gabarito_auto_test = problema_data['gabarito']['resultados_teste_auto']
    
    # Carrega a resolução existente ou inicializa um estado vazio
    resolucao_existente = carregar_resolucao(exercicio_id)
    if resolucao_existente and 'respostas_aluno' in resolucao_existente:
        respostas_aluno = resolucao_existente['respostas_aluno']
    else:
        respostas_aluno = {}

    st.write("Análise dos axiomas de subespaço vetorial:")

    # Usar st.session_state para manter as respostas na interface
    # Inicializa o estado se for a primeira vez ou se o exercício mudou
    if 'current_problem' not in st.session_state or st.session_state.current_problem != exercicio_id:
        st.session_state.current_problem = exercicio_id
        st.session_state.eh_subespaco = respostas_aluno.get('eh_subespaco', None)
        st.session_state.vetor_nulo_check = respostas_aluno.get('vetor_nulo', {}).get('valido', None)
        st.session_state.adicao_check = respostas_aluno.get('adicao', {}).get('valido', None)
        st.session_state.multiplicacao_check = respostas_aluno.get('multiplicacao', {}).get('valido', None)
        st.session_state.multiplicacao_contraexemplo_escalar = respostas_aluno.get('multiplicacao', {}).get('contraexemplo', {}).get('escalar', "1")
        st.session_state.multiplicacao_contraexemplo_vetor = str(respostas_aluno.get('multiplicacao', {}).get('contraexemplo', {}).get('vetor_w_valor', [1,0,0]))
        
        st.session_state.test_results = {
            "vetor_nulo": {"status": None, "message": ""},
            "adicao": {"status": None, "message": ""},
            "multiplicacao": {"status": None, "message": ""}
        }

    # --- Funções auxiliares para teste e feedback ---
    def _display_test_result(test_status, test_message):
        if test_status is True:
            st.success(f"✅ Válido: {test_message}")
        elif test_status is False:
            st.error(f"❌ Inválido: {test_message}")
        else:
            st.warning(f"ℹ️ {test_message}")

    def _run_test_and_display(axiom_name, student_valid, test_counterexample_data=None):
        if student_valid is None:
            st.session_state.test_results[axiom_name] = {"status": None, "message": "Responda 'Sim' ou 'Não' para testar."}
            return

        expected_valid = gabarito_auto_test.get(axiom_name)
        if expected_valid is None:
            st.session_state.test_results[axiom_name] = {"status": None, "message": "Teste automático não configurado para este axioma."}
            return
        
        test_result, message = solver_utils.perform_subspace_axiom_test(
            axiom_name, dados_exercicio, counterexample_data=test_counterexample_data
        )
        
        st.session_state.test_results[axiom_name] = {"status": test_result, "message": message}
        
        # Feedback sobre a resposta do aluno
        if student_valid == expected_valid:
            st.success(f"Sua resposta ('{'Sim' if student_valid else 'Não'}') está correta.")
        else:
            st.error(f"Sua resposta ('{'Sim' if student_valid else 'Não'}') está incorreta. O correto é '{'Sim' if expected_valid else 'Não'}'.")

        _display_test_result(test_result, message)


    # Resposta final do aluno
    st.radio(
        "**Conclusão: O conjunto W é um subespaço de V?**",
        (True, False),
        key='eh_subespaco',
        format_func=lambda x: "Sim" if x else "Não",
        horizontal=True
    )
    st.markdown("---")

    # --- Verificação dos Axiomas ---
    cols = st.columns(3)
    
    with cols[0]:
        st.write("**1. Contém o vetor nulo?**")
        st.radio("Resposta:", (True, False), key='vetor_nulo_check', format_func=lambda x: "Sim" if x else "Não", horizontal=True)
        if st.session_state.vetor_nulo_check is not None:
            _run_test_and_display("vetor_nulo", st.session_state.vetor_nulo_check)

    with cols[1]:
        st.write("**2. É fechado sob adição?**")
        st.radio("Resposta:", (True, False), key='adicao_check', format_func=lambda x: "Sim" if x else "Não", horizontal=True)
        if st.session_state.adicao_check is not None:
            _run_test_and_display("adicao", st.session_state.adicao_check)

    with cols[2]:
        st.write("**3. É fechado sob multiplicação por escalar?**")
        st.radio("Resposta:", (True, False), key='multiplicacao_check', format_func=lambda x: "Sim" if x else "Não", horizontal=True)
        
        # Lógica para contraexemplo de multiplicação se a resposta for "Não" e o gabarito também indicar "Não"
        if st.session_state.multiplicacao_check is False and gabarito_auto_test.get("multiplicacao") is False:
            st.markdown("**Forneça um contraexemplo:**")
            col_escalar, col_vetor = st.columns([1,2])
            with col_escalar:
                st.text_input("Escalar $\lambda$", value=st.session_state.multiplicacao_contraexemplo_escalar, key='multiplicacao_contraexemplo_escalar')
            with col_vetor:
                st.text_input("Vetor $w$ (e.g., [1,0,0])", value=st.session_state.multiplicacao_contraexemplo_vetor, key='multiplicacao_contraexemplo_vetor')
            
            try:
                ce_vetor_str = st.session_state.multiplicacao_contraexemplo_vetor
                # Remove colchetes e divide por vírgula
                ce_vetor_list = [sympy.sympify(x.strip()) for x in ce_vetor_str.strip('[]').split(',') if x.strip()]
                
                counterexample_data = {
                    "escalar": st.session_state.multiplicacao_contraexemplo_escalar,
                    "vetor_w_valor": ce_vetor_list
                }
            except Exception as e:
                st.error(f"Erro ao parsear o vetor de contraexemplo: {e}")
                counterexample_data = None
        else:
            counterexample_data = None

        if st.session_state.multiplicacao_check is not None:
            _run_test_and_display("multiplicacao", st.session_state.multiplicacao_check, counterexample_data=counterexample_data)

    st.markdown("---")
    
    if st.button("Salvar Resolução"):
        resolucao_data = {
            "exercicio_id": exercicio_id,
            "respostas_aluno": {
                "eh_subespaco": st.session_state.eh_subespaco,
                "vetor_nulo": {
                    "valido": st.session_state.vetor_nulo_check,
                    "teste_auto": st.session_state.test_results["vetor_nulo"]
                },
                "adicao": {
                    "valido": st.session_state.adicao_check,
                    "teste_auto": st.session_state.test_results["adicao"]
                },
                "multiplicacao": {
                    "valido": st.session_state.multiplicacao_check,
                    "teste_auto": st.session_state.test_results["multiplicacao"]
                }
            }
        }
        # Adicionar o contraexemplo se foi fornecido e relevante
        if st.session_state.multiplicacao_check is False and gabarito_auto_test.get("multiplicacao") is False:
             resolucao_data["respostas_aluno"]["multiplicacao"]["contraexemplo"] = {
                "escalar": st.session_state.multiplicacao_contraexemplo_escalar,
                "vetor_w_valor": st.session_state.multiplicacao_contraexemplo_vetor # Salva como string para consistência
            }
        
        salvar_resolucao(exercicio_id, resolucao_data)
        st.success(f"Resolução para o exercício '{exercicio_id}' salva com sucesso!")
        st.toast("Não se esqueça de fazer o commit e push no GitHub!")


def renderizar_sistema_homogeneo(problema_data, problem_id="custom"):
    """Renderiza a interface do solver de sistemas lineares."""
    exercicio_id = problem_id
    st.info("O objetivo aqui é escalonar a matriz para encontrar a base do espaço-solução (núcleo).")

    # Inicialização ou recuperação do estado do solver
    if 'solver_state' not in st.session_state or st.session_state.solver_state.get('exercise_id') != exercicio_id:
        equacoes = problema_data['dados_problema']['equacoes']
        matrix, variables, error = parser.parse_system_input(equacoes)
        if error:
            st.error(f"Erro ao analisar o sistema: {error}")
            return
        
        st.session_state.solver_state = {
            'exercise_id': exercicio_id,
            'current_matrix': matrix,
            'variables': variables,
            'matrix_history': [matrix],
            'op_history': [],
            'redo_history': [], # Novo
            'redo_matrix_history': [], # Novo
            'run_id': 0 # Para forçar a atualização dos widgets
        }
    
    solver_state = st.session_state.solver_state
    
    # --- Interface do Solver ---
    st.subheader("Matriz Aumentada")
    matrix = solver_state['current_matrix']
    variables = solver_state['variables']
    num_rows = matrix.shape[0]

    cols = st.columns([0.5] + [1] * matrix.cols)
    for i, col_name in enumerate(["L"] + [f"${v}$" for v in variables] + ["b"]):
        cols[i].markdown(f"**{col_name}**")

    for r in range(num_rows):
        cols = st.columns([0.5] + [1] * matrix.cols)
        cols[0].markdown(f"**$L_{{{r+1}}}$**")
        for c in range(matrix.cols):
            cols[c+1].latex(latex(matrix[r, c]))

    st.markdown("---")

    # --- Controles de Operação ---
    with st.expander("Aplicar Operação de Linha", expanded=True):
        col_op_type, col_params = st.columns([1,2])
        op_type = col_op_type.selectbox("Operação", ["Troca de Linhas", "Multiplicação por Escalar", "Combinação Linear"], key=f"op_select_{solver_state['run_id']}")

        kwargs = {}
        if op_type == "Troca de Linhas":
            r1 = col_params.number_input("Linha A", 1, num_rows, 1, key=f"r1_swap_{solver_state['run_id']}")
            r2 = col_params.number_input("Linha B", 1, num_rows, 2, key=f"r2_swap_{solver_state['run_id']}")
            kwargs = {'r1': r1, 'r2': r2}
        elif op_type == "Multiplicação por Escalar":
            r1 = col_params.number_input("Linha", 1, num_rows, 1, key=f"r_mult_{solver_state['run_id']}")
            alpha = col_params.text_input("Escalar α", "1", key=f"alpha_mult_{solver_state['run_id']}")
            kwargs = {'r': r1, 'alpha': alpha}
        else: # Combinação Linear
            r1 = col_params.number_input("Linha Alvo", 1, num_rows, 1, key=f"r1_comb_{solver_state['run_id']}")
            r2 = col_params.number_input("Linha Base", 1, num_rows, 2, key=f"r2_comb_{solver_state['run_id']}")
            lambda_val = col_params.text_input("Fator λ", "1", key=f"lambda_comb_{solver_state['run_id']}")
            kwargs = {'ri': r1, 'rj': r2, 'factor_lambda': lambda_val}

        if st.button("Aplicar Operação Manual", key=f"apply_manual_op_{solver_state['run_id']}"):
            apply_op(op_type, source="manual", **kwargs)
    
    # --- Navegação e Automação ---
    st.markdown("---")
    st.subheader("Navegação e Automação")
    
    col_nav1, col_nav2, col_auto1, col_auto2, col_auto3 = st.columns(5)
    with col_nav1:
        st.button("Desfazer ↩️", on_click=undo_op, disabled=len(solver_state['op_history']) == 0, key=f"undo_btn_{solver_state['run_id']}")
    with col_nav2:
        st.button("Refazer ↪️", on_click=redo_op, disabled=len(solver_state['redo_history']) == 0, key=f"redo_btn_{solver_state['run_id']}")
    with col_auto1:
        st.button("Próximo Passo", on_click=apply_next_auto_step, key=f"next_step_btn_{solver_state['run_id']}")
    with col_auto2:
        st.button("Resolver para REF", on_click=solve_to_ref_auto, key=f"solve_ref_btn_{solver_state['run_id']}")
    with col_auto3:
        st.button("Resolver para RREF", on_click=solve_all_steps_auto, key=f"solve_rref_btn_{solver_state['run_id']}")

    # --- Histórico de Operações ---
    st.markdown("---")
    if solver_state['op_history']:
        st.write("**Histórico de Operações:**")
        for op in reversed(solver_state['op_history']):
            st.latex(op)

    st.markdown("---")
    if st.button("Salvar Resolução do Sistema", key=f"save_res_btn_{solver_state['run_id']}"):
        final_matrix_str = [str(e) for e in solver_state['current_matrix'].tolist()]
        resolucao_data = {
            "exercicio_id": exercicio_id,
            "respostas_aluno": {
                "tipo": "SISTEMA_HOMOGENEO",
                "historico_operacoes": solver_state['op_history'],
                "matriz_final": final_matrix_str,
                "formato_matriz": solver_state['current_matrix'].shape,
                "variables": solver_state['variables'] # Salvar variáveis para análise futura
            }
        }
        salvar_resolucao(exercicio_id, resolucao_data)
        st.success(f"Resolução para o sistema '{exercicio_id}' salva com sucesso!")
        st.toast("Não se esqueça de fazer o commit e push no GitHub!")

    # --- Ver Resolução Completa (Log) ---
    st.markdown("---")
    with st.expander("Ver Resolução Completa (para copiar)"):
        log_items = solver_utils.generate_solution_log_items(
            solver_state['matrix_history'],
            solver_state['op_history'],
            problema_data['dados_problema']['equacoes'], # User input
            solver_state['variables']
        )
        for item_type, content in log_items:
            if item_type and content:
                if item_type == 'header': st.markdown(content)
                elif item_type == 'matrix_str': st.code(content, language=None)


def render_aluno_view(problem_type, problem_id, problema_data):
    """Chama a função de renderização apropriada com base no tipo de problema."""
    if problem_type == "VERIFICAR_SUBESPACO":
        renderizar_verificacao_subespaco(problema_data, problem_id)
    elif problem_type == "SISTEMA_HOMOGENEO":
        renderizar_sistema_homogeneo(problema_data, problem_id)
    else:
        st.warning(f"O tipo de problema '{problem_type}' ainda não é suportado.")


# --- SIDEBAR ---
st.sidebar.title("Configurações")
st.sidebar.markdown("---")


# --- ÁREA PRINCIPAL ---

st.title("Assistente de Álgebra Linear Interativo")
st.markdown("Bem-vindo ao seu ambiente de trabalho de Álgebra Linear! Aqui você pode explorar, resolver e analisar problemas passo a passo.")

# Seção de entrada de problema customizado
st.subheader("Resolver Problema Customizado")
problem_type_selected = st.selectbox(
    "Selecione o Tipo de Problema:",
    options=["Sistema Linear", "Verificar Subespaço"], # Por enquanto, apenas sistema linear
    key="custom_problem_type_selector"
)

# Inicializa o input do usuário na session_state se não existir
if 'custom_user_input' not in st.session_state:
    st.session_state.custom_user_input = "x+y+z=1\nx-y+z=0\n2x+y-z=2" # Exemplo padrão

if problem_type_selected == "Sistema Linear":
    st.markdown("Cole seu sistema de equações lineares abaixo (uma equação por linha):")
    custom_equations_input = st.text_area("Equações:", value=st.session_state.custom_user_input, height=150, key="custom_system_equations_input")
    
    if st.button("Analisar Sistema Customizado", type="primary"):
        st.session_state.custom_user_input = custom_equations_input # Salva o último input
        matrix, variables, error = parser.parse_system_input(custom_equations_input)
        if error:
            st.error(f"Erro ao analisar o sistema customizado: {error}")
            st.session_state.custom_problem_data = None
        else:
            st.success("Sistema customizado processado!")
            st.info(f"Variáveis detectadas: {', '.join(variables)}")
            # Estrutura como um 'exercicio' para ser compatível com as funções de renderização
            st.session_state.custom_problem_data = {
                "id": "custom_linear_system",
                "titulo": "Sistema Linear Customizado",
                "enunciado": f"Sistema inserido:\n```\n{custom_equations_input}\n```",
                "tipo": "SISTEMA_HOMOGENEO", # Usar o renderer existente
                "dados_problema": { # Mudança para 'dados_problema' para diferenciar de 'dados_exercicio' do JSON
                    "equacoes": custom_equations_input
                },
                "gabarito": { # Gabarito genérico para custom, pode ser a análise final
                    "eh_subespaco": False # Não aplicável aqui
                }
            }
        st.rerun() # Forçar rerun para exibir o solver
        
elif problem_type_selected == "Verificar Subespaço":
    st.markdown("Cole a definição do subespaço $W$ e do espaço $V$ (formato específico):")
    st.warning("Este módulo ainda não está implementado para problemas customizados.")


# Exibir o solver interativo se houver dados de problema customizado
if 'custom_problem_data' in st.session_state and st.session_state.custom_problem_data is not None:
    problem_type = st.session_state.custom_problem_data['tipo']
    problem_id = st.session_state.custom_problem_data['id']
    problema_data = st.session_state.custom_problem_data
    
    st.markdown("---")
    st.subheader(f"Análise do Problema Customizado: {problema_data['titulo']}")
    st.markdown(problema_data['enunciado'], unsafe_allow_html=True)
    st.markdown("---")

    render_aluno_view(problem_type, problem_id, problema_data)

    # A área de gabarito permanece para o aluno consultar, se desejar.
    with st.expander("Ver Gabarito (para referência, se aplicável)"):
        st.write("Não há gabarito pré-definido para problemas customizados. A análise automática é exibida no solver.")

else:
    st.info("Selecione um tipo de problema e cole sua fórmula para começar a análise.")
