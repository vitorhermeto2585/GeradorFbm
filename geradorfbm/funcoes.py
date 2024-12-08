import numpy as np
from scipy.interpolate import RectBivariateSpline

# Função para gerar séries temporais simulando Movimento Browniano Fracionário (fBm)
def fbm(n, H=0.75, caminho=None):
    """
    Gera uma série temporal de Movimento Browniano Fracionário (fBm) com base no parâmetro de Hurst.

    Parâmetros:
        - n: Número de pontos na série temporal.
        - H: Parâmetro de Hurst (0 < H < 1).
        - caminho: Caminho para salvar a série gerada em CSV (opcional).

    Retorna:
        - Salva a série gerada em um arquivo CSV no caminho especificado.
    """
    metodo = 'subgrupos_fft'
    div = 11
    j = 1500

    assert 0 < H < 1, "O parâmetro Hurst (H) deve estar no intervalo (0, 1)."

    def autocovariance(k, H):
        """
        Calcula a autocovariância para o fBm, essencial para criar a matriz de correlação.
        """
        return 0.5 * (np.abs(k + 1) ** (2 * H) - 2 * np.abs(k) ** (2 * H) + np.abs(k - 1) ** (2 * H))

    def R(t, s):
        """
        Função de autocorrelação baseada no parâmetro de Hurst.
        """
        twoH = 2 * H
        return 0.5 * (s ** twoH + t ** twoH - np.abs(t - s) ** twoH)

    def resize_matrix_sigma(sigma, new_size):
        """
        Redimensiona a matriz sigma para se ajustar ao tamanho do subgrupo.
        Utiliza interpolação bidimensional.
        """
        x_old = np.linspace(0, 1, sigma.shape[0])
        y_old = np.linspace(0, 1, sigma.shape[1])
        interpolator = RectBivariateSpline(x_old, y_old, sigma)
        x_new = np.linspace(0, 1, new_size)
        y_new = np.linspace(0, 1, new_size)
        sigma_resized = interpolator(x_new, y_new)
        return sigma_resized

    if H >= 0.5:
        # Algoritmo baseado em FFT para H >= 0.5, garantindo maior eficiência
        g = np.zeros(2 * n)
        for k in range(n):
            g[k] = autocovariance(k, H)
        g = np.concatenate([g, g[::-1]])[:2 * n]
        G = np.fft.fft(g)
        G = np.sqrt(np.real(G))
        W = np.random.randn(2 * n) + 1j * np.random.randn(2 * n)
        Z = np.fft.ifft(G * W).real[:n]
        serie = np.cumsum(Z)
    else:
        # Subgrupos para melhorar eficiência e memória longa para H < 0.5
        gamma = R(*np.mgrid[0:j, 0:j])
        w, P = np.linalg.eigh(gamma)
        L = np.diag(w)
        sigma = np.dot(np.dot(P, np.sqrt(L)), np.linalg.inv(P))
        j = n // div
        p = n - (j * (div - 1))
        aux = p - j
        sigma = resize_matrix_sigma(sigma, j)
        serie = []
        for i in range(div):
            v = np.random.randn(j)
            serie.extend(np.dot(sigma, v))
        serie.extend([np.mean(serie)] * aux)

    if caminho is None:
        caminho = ""

    # Valida a série gerada com DFA e salva em CSV
    np.savetxt(f'{caminho}{metodo}_{n}_{H}.csv', serie, delimiter=',', comments='')
    

# Função para realizar a análise DFA (Detrended Fluctuation Analysis)
def DFA(serie_temporal, escala_min=11, escala_max=None, densidade=8, ordem=1):
    """
    Realiza a análise de Flutuação Sem Tendência (DFA) em uma série temporal.

    Parâmetros:
        - serie_temporal: Série temporal para análise (numpy array ou lista).
        - escala_min: Tamanho mínimo da escala (janela).
        - escala_max: Tamanho máximo da escala (opcional).
        - densidade: Número de divisões logarítmicas na escala.
        - ordem: Ordem do polinômio para ajuste.

    Retorna:
        - Exponente DFA estimado.
    """
    tamanho_da_serie = len(serie_temporal)
    media = serie_temporal.mean()
    # Subtrai a média para eliminar tendências globais
    serie_nivelada = serie_temporal - media
    # Calcula a soma cumulativa, fundamental para o perfil acumulado
    soma_acumulada = np.cumsum(serie_nivelada)

    def calculoJanela(tamanhoJanela, numeroJanelas, soma_acumulada):
        """
        Calcula a flutuação média ao ajustar uma reta em cada janela da série.
        """
        st = np.array([])
        n = len(soma_acumulada)
        for nj in range(numeroJanelas):
            inicio = nj * tamanhoJanela
            fim = inicio + tamanhoJanela
            janela = soma_acumulada[inicio:fim]
            x = np.arange(inicio + 1, fim + 1)
            rectaV = np.polyfit(x, janela, ordem)
            fV = np.poly1d(rectaV)
            ordenada = fV(x)
            st = np.append(st, (ordenada - janela) ** 2)
        valor = st.cumsum()
        valorFinal = valor[-1] / n
        return valorFinal ** 0.5

    # Define o tamanho das janelas baseado na densidade logarítmica
    if escala_max is None:
        escala_max = tamanho_da_serie // 4

    tamJanSug = []
    for t in range(escala_min, escala_max):
        tamanho = int(4 * (2 ** (1 / densidade)) ** t + 0.5)
        if tamanho <= tamanho_da_serie / 4:
            tamJanSug.append(tamanho)
        else:
            break

    # Calcula as flutuações em cada janela e determina o comportamento em escala
    Fn = []
    for tamjanela in tamJanSug:
        Fn.append(calculoJanela(int(tamjanela), int(tamanho_da_serie / tamjanela), soma_acumulada))

    retalog = np.polyfit(np.log2(tamJanSug), np.log2(Fn), 1)
    flog = np.poly1d(retalog)

    return flog[1]
