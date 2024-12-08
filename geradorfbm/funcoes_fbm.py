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
