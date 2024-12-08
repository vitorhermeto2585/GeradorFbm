# GeradorFbm

Uma biblioteca Python desenvolvida para a geração de séries temporais utilizando Movimento Browniano Fracionário (fBm) e análise DFA (Detrended Fluctuation Analysis). A biblioteca oferece métodos precisos e eficientes para simular e validar propriedades estatísticas de séries temporais com diferentes características de memória.

---

## Visão Geral

### Objetivos e Funcionalidades
O **GeradorFbm** combina métodos avançados de geração de séries temporais e análise estatística para:
- **Simulação precisa de séries com memória longa e curta**:
  - Memória longa (H > 0.5): Séries com persistência nas flutuações.
  - Memória curta ou anti-persistente (H < 0.5): Séries com reversão frequente nas flutuações.
- **Validação robusta das propriedades das séries geradas**:
  - Implementação do método DFA para calcular o expoente de Hurst, um indicador crucial da memória das séries temporais.

### Por que usar esta biblioteca?
- **Flexibilidade**: Permite configurar o parâmetro de Hurst para gerar séries com características de memória ajustáveis.
- **Eficiência Computacional**: Utiliza métodos otimizados, como subgrupos e transformadas rápidas de Fourier (FFT), garantindo desempenho para séries grandes.
- **Precisão Estatística**: Preserva as propriedades desejadas das séries geradas, alinhadas ao expoente de Hurst configurado.


---

### Principais Funcionalidades

1. **`fbm_final`**: Gera uma série temporal simulando um Movimento Browniano Fracionário com base no parâmetro de Hurst.
   - **Parâmetros**:
     - `n` (int): Número de pontos na série temporal.
     - `H` (float): Parâmetro de Hurst (0 < H < 1).
     - `caminho` (str): Caminho para salvar a série gerada em CSV (opcional).
   - **Retorno**: Salva a série gerada em um arquivo CSV

2. **`DFA`**: Analisa séries temporais para calcular o expoente de Hurst.
   - **Parâmetros**:
     - `serie_temporal` (list ou numpy array): Série temporal para análise.
     - `escala_min` (int): Tamanho mínimo das janelas (default: 11).
     - `escala_max` (int): Tamanho máximo das janelas (default: None).
     - `densidade` (int): Número de divisões logarítmicas na escala (default: 8).
     - `ordem` (int): Ordem do polinômio usado para ajuste (default: 1).
   - **Retorno**: O expoente DFA (float), que estima a memória longa da série.


---

## Instalação

Para instalar a biblioteca diretamente do GitHub, execute:

```bash
pip install git+https://github.com/vitorhermeto2585/GeradorFbm.git
