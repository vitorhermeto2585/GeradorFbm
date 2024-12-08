import numpy as np

# Função para realizar a análise DFA (Detrended Fluctuation Analysis)
def dfa(serie_temporal, escala_min=11, escala_max=None, densidade=8, ordem=1):
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
