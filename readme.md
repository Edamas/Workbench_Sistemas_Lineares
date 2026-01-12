# Workbench de Sistemas Lineares
**Álgebra Linear — Projeto Interativo**

---

## 1. Visão Geral

Este projeto é um **ambiente interativo de aprendizado** (workbench), construído em Streamlit, projetado para que estudantes possam resolver sistemas de equações lineares através do **Método de Eliminação de Gauss e Gauss-Jordan**. A ferramenta foca em um processo de resolução passo a passo, manual e consciente, onde cada operação elementar é registrada e seus efeitos são visualizados em tempo real.

O objetivo é **educar e capacitar**, não apenas automatizar o cálculo. O estudante mantém controle total do processo, enquanto o sistema oferece validação matemática e feedback instantâneo.

---

## 2. Funcionalidades Principais

### 2.1 Resolução Interativa
O coração da ferramenta é a tabela de escalonamento interativa, que permite ao usuário:
- **Digitar ou gerar** um sistema de equações lineares.
- **Visualizar** a matriz aumentada e as equações correspondentes em tempo real.
- **Aplicar manualmente** as três operações elementares de linha (Troca, Escala, Combinação).
- **Receber feedback instantâneo** sobre a validade de cada operação e o estado da matriz (pivôs, linhas nulas, inconsistências).
- **Desfazer e refazer** operações, navegando livremente pelo histórico de resolução.

### 2.2 Gerador de Sistemas Aleatórios
Para facilitar a prática e o teste de diferentes cenários, a ferramenta inclui um gerador de sistemas configurável com as seguintes opções:
- **Número de Variáveis:** De 1 a 10.
- **Notação das Variáveis:** `x, y, z...` ou `x₁, x₂, x₃...`.
- **Tipos de Coeficientes:** Inclusão opcional de frações, números complexos e raízes.
- **Tipo de Solução:** Garante que o sistema gerado seja Determinado (SPD), Indeterminado (SPI) ou Impossível (SI).
- **Formato da Multiplicação:** Controle sobre a exibição da multiplicação (implícita, `*` ou `⋅`) e o uso de parênteses.

### 2.3 Análise e Relatório
Ao final do processo de escalonamento, o workbench oferece:
- **Classificação automática** do sistema (SPD, SPI, SI) com base no Teorema de Rouché-Capelli.
- **Exibição da solução final** do sistema.
- **Um relatório completo e detalhado** ("Ver Resolução Completa") que pode ser copiado, mostrando o sistema original, a matriz inicial e todas as etapas da resolução.

---

## 3. Roadmap de Melhorias Propostas

As funcionalidades abaixo foram planejadas para enriquecer a experiência pedagógica e de usabilidade da ferramenta.

#### 1. Badge Fixo de Estado da Matriz
Implementar um badge visualmente proeminente e fixo na interface que indique o estado atual da matriz em tempo real:
- ❌ **Não escalonada**
- ⚠️ **Escalonada**
- ✅ **Escalonada Reduzida**
Lembrando que, para o material de ensino ministrado pelo Professor do IME:
Um sistema linear (S) é escalonado se:
1. a primeira variável presente em uma linha estiver a direita da primeira variável da linha
superior;
2. a primeira variável presente na equação de uma linha não estiver presente nas linhas acima;
3. a primeira variável de cada equação é seguida do coeficiente 1;
4. qualquer linha nula estiver abaixo de qualquer linha não nula.

#### 2. Destaque Visual do Pivô
Para ancorar o conceito de pivô, serão adicionadas melhorias visuais na tabela interativa:
- **Fundo Sombreado:** A coluna do pivô atual terá um leve sombreamento.
- **Ícone de Alvo (🎯):** A linha do pivô será marcada com um ícone para fácil identificação.
- **Tooltip Explicativo:** Ao passar o mouse sobre o pivô, um tooltip informará: “Este é o pivô: o primeiro coeficiente não nulo da linha.”

#### 3. Padronização da Linguagem Técnica
A terminologia na interface será refinada para alinhar-se com a literatura clássica de Álgebra Linear:
- **Título da Seção:** "Operações Elementares de Linha".
- **Nomes das Operações:** "Troca", "Escala" e "Combinação".

#### 4. Melhoria na Exibição LaTeX
O histórico de operações na barra lateral será aprimorado para renderizar uma saída LaTeX pura e sem ruídos, padronizando frações no formato `\frac{a}{b}`.

#### 5. Modo "Professor" vs. "Aluno"
Uma nova opção de alternância será criada para adaptar a ferramenta a diferentes contextos de uso:
- **Modo Aluno (Padrão):** Exibe todas as observações automáticas, alertas e dicas (ex: “pivô fora de ordem”, “linha nula encontrada”).
- **Modo Professor:** Oculta as observações automáticas, permitindo que a ferramenta seja usada para avaliação, onde apenas o resultado final da validação é exibido.

#### 6. Coluna "Justificativa Matemática"
Uma coluna opcional será adicionada à tabela de escalonamento para que o usuário possa descrever o **objetivo matemático** de cada operação. Exemplo: "Zerar o coeficiente da variável x₁ na linha 3".

#### 7. Relatório Final Exportável
A seção "Ver Resolução Completa" será expandida com botões para exportar o relatório final nos seguintes formatos:
- 📄 **Exportar Markdown**
- 📄 **Exportar LaTeX**
- 📄 **Exportar PDF** (gerado a partir do LaTeX)

#### 8. Nomear o Método
O nome do método utilizado será explicitamente exibido na interface principal (ex: "Método de Eliminação de Gauss / Gauss-Jordan") para reforçar o conceito teórico.

---

## 4. Itens que Precisam de Atenção

- A implementação do **Roadmap de Melhorias** representa um esforço de desenvolvimento significativo. Cada item será implementado de forma iterativa.
- A consistência da renderização **LaTeX** em toda a aplicação, especialmente no histórico de operações, precisa ser revisada para garantir que apenas o formato matemático seja exibido, sem texto ou formatação adicional indesejada.

---

## 6. Diretrizes Operacionais (Backups)

Para garantir uma gestão eficiente de versões e evitar o excesso de arquivos, as seguintes diretrizes para backups serão estritamente seguidas:

-   **Backup Inicial por Solicitação:** Um backup dos arquivos relevantes será criado **apenas uma vez** no início de cada nova solicitação do usuário que envolva modificações no código. Este backup representa o estado original dos arquivos antes de qualquer alteração para aquela solicitação específica.
-   **Nomenclatura Padrão:** Os arquivos de backup serão salvos no diretório `backups/` e seguirão o formato `[nome_arquivo_original]_[carimbo_de_data_e_hora].py` (ex: `app.py_YYYY-MM-DD-HH-MM.py`).
-   **Sem Backups Intermediários:** Durante a execução de uma única solicitação do usuário, **não serão criados backups adicionais**. Todas as modificações serão aplicadas diretamente ao arquivo de trabalho.
-   **Ponto de Reversão:** O backup inicial de uma solicitação serve como um ponto de reversão, caso o usuário não aprove as alterações ou deseje restaurar o estado anterior.
-   **Retenção de Backups:** Backups de solicitações anteriores serão gerenciados para reter apenas as versões mais significativas, conforme determinado pela aprovação explícita do usuário de grandes mudanças ou entrega final. Backups intermediários de uma mesma solicitação (se criados por engano) serão removidos para reduzir o volume.
-   **Foco no Arquivo Ativo:** Embora erros em arquivos de backup sejam notados, o foco principal de correção e desenvolvimento será sempre o arquivo de trabalho ativo (`app.py`, `parser.py`, etc.).
-   **Objetividade e Eficiência:** O trabalho será conduzido com o máximo de objetividade, praticidade, rapidez e evitando loops desnecessários.

---

## 5. Princípios Norteadores
- **Clareza Matemática > Automação:** O foco é no entendimento do processo.
- **Processo > Resposta:** O caminho da resolução é mais importante que a solução final.
- **Aprendizado Ativo:** O usuário está no controle de cada decisão.
- **Rigor sem Perder Usabilidade:** A interface deve ser intuitiva, mas fiel à terminologia e aos processos formais.