
# Análise de Cariótipos

Este projeto realiza uma análise básica de cariótipos utilizando as bibliotecas BioPython e matplotlib. Ele foi desenvolvido por Mariana Silva de Oliveira.

## Descrição

O script `analise_cariotipos.py` lê sequências de DNA de um arquivo FASTA, calcula o tamanho e o conteúdo GC de cada cromossomo, e gera um gráfico com essas informações.

## Requisitos

- Python 3.6+
- BioPython
- matplotlib

## Instalação

1. Clone este repositório:
   ```
   git clone https://github.com/ezrafchev/analise-de-cariotipo-by-mariana.git
   cd analise-de-cariotipo-by-mariana
   ```

2. Instale as dependências:
   ```
   pip install biopython matplotlib
   ```

## Uso

1. Prepare seu arquivo FASTA com as sequências dos cromossomos. Um exemplo está incluído no arquivo `exemplo_cariotipos.fasta`.

2. Execute o script:
   ```
   python3 analise_cariotipos.py
   ```

3. O script gerará um arquivo `resultados_cariotipos.png` com os gráficos de análise.

## Estrutura do Projeto

- `analise_cariotipos.py`: Script principal para análise de cariótipos.
- `exemplo_cariotipos.fasta`: Arquivo FASTA de exemplo com sequências de cromossomos.
- `resultados_cariotipos.png`: Gráfico gerado com os resultados da análise.

## Interpretação dos Resultados

O gráfico gerado contém duas partes:

1. Tamanho dos cromossomos: Mostra o comprimento de cada cromossomo em pares de base.
2. Conteúdo GC dos cromossomos: Exibe a porcentagem de bases G e C em cada cromossomo.

Estas informações são úteis para comparar características básicas entre diferentes cromossomos ou entre espécies.

## Limitações e Melhorias Futuras

- O script atual realiza apenas análises básicas. Futuras versões poderiam incluir:
  - Identificação de regiões repetitivas
  - Análise de genes e elementos regulatórios
  - Comparações entre espécies
- A visualização poderia ser melhorada com gráficos interativos

## Contribuições

Contribuições são bem-vindas! Por favor, abra uma issue para discutir mudanças maiores antes de enviar um pull request.

## Créditos

Este projeto foi desenvolvido por Mariana Silva de Oliveira.

## Licença

Este projeto está sob a licença MIT. Veja o arquivo `LICENSE` para mais detalhes.
