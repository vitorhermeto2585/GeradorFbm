# __init__.py

# Define a versão do pacote
__version__ = "0.1.0"

# Importa as funções principais do pacote
from .funcoes import dfa, fbm

# Lista de elementos disponíveis ao importar o pacote
__all__ = ["dfa", "fbm"]
