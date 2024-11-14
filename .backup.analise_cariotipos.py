
import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt

def analisar_cariotipos(arquivo_fasta):
    sequencias = list(SeqIO.parse(arquivo_fasta, "fasta"))
    resultados = []
    
    for seq in sequencias:
        tamanho = len(seq)
        gc_content = gc_fraction(seq.seq) * 100  # Convertendo para porcentagem
        resultados.append({
            'nome': seq.id,
            'tamanho': tamanho,
            'gc_content': gc_content
        })
    
    return resultados

def visualizar_resultados(resultados):
    nomes = [r['nome'] for r in resultados]
    tamanhos = [r['tamanho'] for r in resultados]
    gc_contents = [r['gc_content'] for r in resultados]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    
    ax1.bar(nomes, tamanhos)
    ax1.set_title('Tamanho dos cromossomos')
    ax1.set_xlabel('Cromossomos')
    ax1.set_ylabel('Tamanho (pb)')
    ax1.tick_params(axis='x', rotation=45)
    
    ax2.bar(nomes, gc_contents)
    ax2.set_title('Conteúdo GC dos cromossomos')
    ax2.set_xlabel('Cromossomos')
    ax2.set_ylabel('Conteúdo GC (%)')
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig('resultados_cariotipos.png')
    plt.close()

if __name__ == "__main__":
    arquivo_fasta = "exemplo_cariotipos.fasta"
    if not os.path.exists(arquivo_fasta):
        print(f"Erro: O arquivo {arquivo_fasta} não foi encontrado.")
    else:
        resultados = analisar_cariotipos(arquivo_fasta)
        visualizar_resultados(resultados)
        print("Análise concluída. Resultados salvos em 'resultados_cariotipos.png'")
