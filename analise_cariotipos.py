
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt

def analisar_cariotipos(arquivo_fasta):
    # Carregar sequências do arquivo FASTA
    sequencias = list(SeqIO.parse(arquivo_fasta, "fasta"))
    
    # Análise básica
    resultados = []
    for seq_record in sequencias:
        tamanho = len(seq_record.seq)
        gc_content = (seq_record.seq.count('G') + seq_record.seq.count('C')) / tamanho * 100
        resultados.append({
            'nome': seq_record.id,
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
    arquivo_fasta = "exemplo_cariotipos.fasta"  # Substitua pelo nome do seu arquivo FASTA
    resultados = analisar_cariotipos(arquivo_fasta)
    visualizar_resultados(resultados)
    print("Análise concluída. Resultados salvos em 'resultados_cariotipos.png'")
