
import cv2
import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
from scipy import ndimage

def carregar_imagem(caminho):
    return cv2.imread(caminho, cv2.IMREAD_GRAYSCALE)

def pre_processar_imagem(imagem):
    # Aplicar filtro gaussiano para reduzir ruído
    imagem_suavizada = cv2.GaussianBlur(imagem, (5, 5), 0)
    
    # Aplicar limiarização adaptativa
    imagem_binaria = cv2.adaptiveThreshold(imagem_suavizada, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY_INV, 11, 2)
    
    # Aplicar operações morfológicas para melhorar a segmentação
    kernel = np.ones((3,3), np.uint8)
    imagem_processada = cv2.morphologyEx(imagem_binaria, cv2.MORPH_CLOSE, kernel, iterations=2)
    
    return imagem_processada

def segmentar_cromossomos(imagem_processada):
    # Rotular componentes conectados
    rotulos = measure.label(imagem_processada)
    
    # Filtrar componentes muito pequenos ou muito grandes
    tamanhos = np.bincount(rotulos.ravel())
    mascara_tamanho = (tamanhos > 100) & (tamanhos < 10000)
    rotulos_filtrados = mascara_tamanho[rotulos]
    
    return rotulos_filtrados

def analisar_cromossomos(rotulos_filtrados):
    propriedades = measure.regionprops(rotulos_filtrados)
    
    resultados = []
    for prop in propriedades:
        # Calcular características do cromossomo
        area = prop.area
        perimetro = prop.perimeter
        comprimento = prop.major_axis_length
        largura = prop.minor_axis_length
        
        # Calcular índice centromérico (simplificado)
        y, x = prop.centroid
        indice_centrometrico = y / prop.bbox[2]  # Posição Y do centróide / altura do cromossomo
        
        resultados.append({
            'area': area,
            'perimetro': perimetro,
            'comprimento': comprimento,
            'largura': largura,
            'indice_centrometrico': indice_centrometrico
        })
    
    return resultados

def visualizar_resultados(imagem_original, rotulos_filtrados, resultados):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
    
    # Imagem original
    ax1.imshow(imagem_original, cmap='gray')
    ax1.set_title('Imagem Original')
    ax1.axis('off')
    
    # Imagem segmentada
    ax2.imshow(rotulos_filtrados, cmap='nipy_spectral')
    ax2.set_title('Cromossomos Segmentados')
    ax2.axis('off')
    
    plt.tight_layout()
    plt.savefig('resultados_segmentacao.png')
    plt.close()
    
    # Gráfico de dispersão: Comprimento vs. Índice Centromérico
    plt.figure(figsize=(10, 6))
    comprimentos = [r['comprimento'] for r in resultados]
    indices_centrometricos = [r['indice_centrometrico'] for r in resultados]
    plt.scatter(comprimentos, indices_centrometricos)
    plt.xlabel('Comprimento')
    plt.ylabel('Índice Centromérico')
    plt.title('Comprimento vs. Índice Centromérico dos Cromossomos')
    plt.savefig('grafico_comprimento_vs_indice.png')
    plt.close()

def main(caminho_imagem):
    imagem_original = carregar_imagem(caminho_imagem)
    imagem_processada = pre_processar_imagem(imagem_original)
    rotulos_filtrados = segmentar_cromossomos(imagem_processada)
    resultados = analisar_cromossomos(rotulos_filtrados)
    visualizar_resultados(imagem_original, rotulos_filtrados, resultados)
    
    print(f"Análise concluída. {len(resultados)} cromossomos identificados.")
    print("Resultados salvos em 'resultados_segmentacao.png' e 'grafico_comprimento_vs_indice.png'")

if __name__ == "__main__":
    caminho_imagem = "exemplo_cariotipo.jpg"  # Substitua pelo caminho da sua imagem
    main(caminho_imagem)
