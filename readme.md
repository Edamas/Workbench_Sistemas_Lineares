# Workbench de Sistemas Lineares

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Streamlit](https://img.shields.io/badge/Streamlit-1.29-FF4B4B.svg)](https://streamlit.io)
[![GitHub issues](https://img.shields.io/github/issues/Edamas/Workbench_Sistemas_Lineares)](https://github.com/Edamas/Workbench_Sistemas_Lineares/issues)

**Um ambiente interativo para o aprendizado prático e visual do escalonamento de matrizes e da resolução de sistemas lineares, inspirado no curso de verão de Álgebra Linear do IME-USP.**

---

### **[Acesse a demonstração ao vivo aqui!]**(link-para-o-streamlit-cloud-aqui)

---

## 1. Visão Geral

Este projeto é um **ambiente de aprendizado interativo** (workbench) construído em Streamlit, projetado para que estudantes de Álgebra Linear possam resolver sistemas de equações através do **Método de Eliminação de Gauss e Gauss-Jordan**. A ferramenta foca em um processo de resolução passo a passo, manual e consciente, onde cada operação elementar é aplicada pelo usuário e seus efeitos são visualizados em tempo real.

O objetivo é **capacitar e educar**, não apenas automatizar o cálculo. O estudante mantém controle total do processo, enquanto o sistema oferece validação matemática, feedback instantâneo e ferramentas para facilitar a prática.

## 2. O Curso de Verão IME-USP (2026)

A inspiração para este projeto nasceu da experiência no curso de verão de **Álgebra Linear** de 2026, ministrado pelo Professor **Heitor [Sobrenome a ser confirmado]** no Instituto de Matemática e Estatística da Universidade de São Paulo (IME-USP).

O curso enfatiza a compreensão profunda dos conceitos fundamentais, como o escalonamento de matrizes e o Teorema de Rouché-Capelli, que são a base para a análise e resolução de sistemas lineares. Esta ferramenta busca digitalizar e enriquecer a experiência de aprendizado proposta em sala de aula, permitindo que alunos e entusiastas possam praticar e visualizar esses conceitos de forma dinâmica.

## 3. Funcionalidades Principais

- **Resolução Interativa:** Aplique as três operações elementares de linha (Troca, Multiplicação por Escalar, Combinação Linear) e veja o impacto imediato na matriz aumentada.
- **Gerador de Sistemas:** Crie sistemas de equações aleatórios, configurando o número de variáveis, o tipo de coeficientes (reais, complexos, frações, raízes) e o tipo de solução (SPD, SPI, SI).
- **Análise em Tempo Real:** Receba feedback instantâneo sobre o estado da matriz (escalonada, escalonada reduzida), a posição dos pivôs e a consistência do sistema.
- **Histórico de Operações:** Desfaça e refaça passos, permitindo a exploração de diferentes estratégias de resolução sem perder o trabalho anterior.
- **Relatório de Resolução:** Ao final, obtenha um relatório completo com todas as etapas, desde o sistema original até a solução, pronto para ser copiado e compartilhado.

## 4. Como Usar

Para executar este projeto localmente, siga os passos abaixo:

1.  **Clone o repositório:**
    ```bash
    git clone https://github.com/Edamas/Workbench_Sistemas_Lineares.git
    cd Workbench_Sistemas_Lineares
    ```

2.  **Crie um ambiente virtual (recomendado):**
    ```bash
    python -m venv .venv
    source .venv/bin/activate  # No Windows, use `.venv\Scripts\activate`
    ```

3.  **Instale as dependências:**
    ```bash
    pip install -r requirements.txt
    ```

4.  **Execute a aplicação Streamlit:**
    ```bash
    streamlit run app.py
    ```

## 5. Próximos Passos e Roadmap

A visão para o futuro deste projeto é torná-lo uma ferramenta pedagógica ainda mais robusta e versátil.

- **Publicação na Streamlit Community Cloud:** O próximo passo imediato é disponibilizar a aplicação online para acesso público.
- **Melhorias na Interface:** Aprimorar a experiência do usuário com destaques visuais para pivôs, uma linguagem técnica mais padronizada e uma melhor renderização de equações em LaTeX.
- **Modos de Uso ("Aluno" e "Professor"):** Implementar um modo "Professor" que oculte as dicas automáticas, permitindo que a ferramenta seja usada para avaliações ou demonstrações.
- **Exportação de Relatórios:** Adicionar a funcionalidade de exportar a resolução completa para formatos como Markdown, LaTeX e PDF.
- **Aplicações Práticas:** Criar uma seção com exemplos de sistemas lineares aplicados a problemas reais em áreas como engenharia, física e ciência da computação.

## 6. Como Contribuir

Este é um projeto de código aberto e acadêmico. Sugestões, ideias e contribuições são muito bem-vindas! Se você tem uma ideia para uma nova funcionalidade, uma aplicação interessante ou encontrou um bug, sinta-se à vontade para:

-   **Abrir uma [Issue](https://github.com/Edamas/Workbench_Sistemas_Lineares/issues)** para discutir sua ideia.
-   **Enviar um [Pull Request](https://github.com/Edamas/Workbench_Sistemas_Lineares/pulls)** com suas melhorias.

## 7. Licença

Este projeto está licenciado sob a **Licença MIT**. Veja o arquivo [LICENSE](LICENSE) para mais detalhes.
