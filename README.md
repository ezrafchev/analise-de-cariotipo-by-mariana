# Análise de Cariótipos

Este projeto realiza uma análise detalhada de cariótipos utilizando as bibliotecas BioPython, NumPy e Matplotlib. Ele foi desenvolvido por Mariana Silva de Oliveira.

## Descrição

O script `analise_cariotipos.py` lê sequências de DNA de um arquivo FASTA e realiza as seguintes análises:

1. Cálculo do tamanho dos cromossomos
2. Análise do conteúdo GC e AT
3. Cálculo do peso molecular
4. Identificação de regiões repetitivas
5. Simulação de densidade de genes
6. Identificação de ilhas CpG
7. Cálculo da complexidade da sequência

O script gera um gráfico com várias visualizações dessas análises.

## Requisitos

- Python 3.6+
- BioPython
- NumPy
- Matplotlib

## Instalação

1. Clone este repositório:
   ```
   git clone https://github.com/ezrafchev/analise-de-cariotipo-by-mariana.git
   cd analise-de-cariotipo-by-mariana
   ```

2. Instale as dependências:
   ```
   pip install biopython numpy matplotlib
   ```

## Uso

1. Prepare seu arquivo FASTA com as sequências dos cromossomos. Um exemplo está incluído no arquivo `exemplo_cariotipos.fasta`.

2. Execute o script:
   ```
   python3 analise_cariotipos.py
   ```

3. O script gerará um arquivo `resultados_cariotipos.png` com os gráficos de análise e imprimirá resultados detalhados no console.

## Estrutura do Projeto

- `analise_cariotipos.py`: Script principal para análise de cariótipos.
- `exemplo_cariotipos.fasta`: Arquivo FASTA de exemplo com sequências de cromossomos.
- `resultados_cariotipos.png`: Gráfico gerado com os resultados da análise.

## Interpretação dos Resultados

O gráfico gerado contém seis visualizações:

1. Tamanho dos cromossomos: Mostra o comprimento de cada cromossomo em pares de base.
2. Composição de bases: Exibe a porcentagem de bases GC e AT em cada cromossomo.
3. Relação entre conteúdo GC e tamanho: Gráfico de dispersão mostrando a relação entre o conteúdo GC e o tamanho dos cromossomos.
4. Densidade de genes (simulada): Mostra a densidade simulada de genes em cada cromossomo.
5. Complexidade da sequência: Exibe a complexidade média da sequência para cada cromossomo.
6. Ideograma do cariótipo: Uma representação visual simplificada do cariótipo.

Além disso, o script imprime informações detalhadas sobre cada cromossomo, incluindo o número de regiões repetitivas e ilhas CpG identificadas.

## Limitações e Melhorias Futuras

- A simulação de densidade de genes é uma aproximação e não reflete dados reais.
- A identificação de regiões repetitivas e ilhas CpG pode ser refinada com algoritmos mais sofisticados.
- Futuros desenvolvimentos poderiam incluir:
  - Análise de elementos regulatórios
  - Comparações entre espécies
  - Integração com bancos de dados genômicos

## Contribuições

Contribuições são bem-vindas! Por favor, abra uma issue para discutir mudanças maiores antes de enviar um pull request.

## Créditos

Este projeto foi desenvolvido por Mariana Silva de Oliveira.

## Licença

Este projeto está sob a licença MIT. Veja o arquivo `LICENSE` para mais detalhes.
