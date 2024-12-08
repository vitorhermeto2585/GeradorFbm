from setuptools import setup, find_packages

setup(
    name="GeradorFbm",  # Nome do pacote
    version="0.1.0",  # Versão inicial
    author="Vitor Hermeto",
    author_email="vitorolhermeto@gmail.com",  # Substitua pelo seu email
    description="Biblioteca para gerar Movimento Browniano Fracionário (fBm)",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/vitorhermeto2585/GeradorFbm",
    packages=find_packages(),  # Detecta automaticamente os pacotes
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",  # Versão mínima do Python
    install_requires=[
        "numpy>=1.18.0",          # Para cálculos matemáticos
        "scipy>=1.4.0",           # Para interpolação e funções científicas
    ],
)
