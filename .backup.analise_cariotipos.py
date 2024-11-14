
import os
from Bio import SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt

def analisar_cariotipos(arquivo_fasta):
    sequencias = list(SeqIO.parse(arquivo_fasta, "fasta"))
    resultados = []
    
    for seq in sequencias:
        tamanho = len(seq)
        gc_content = GC(seq.seq)
        at_content = 100 - gc_content
        resultados.append({
            'nome': seq.id,
            'tamanho': tamanho,
            'gc_content': gc_content,
            'at_content': at_content
        })
    
    return resultados

def visualizar_resultados(resultados):
    nomes = [r['nome'] for r in resultados]
    tamanhos = [r['tamanho'] for r in resultados]
    gc_contents = [r['gc_content'] for r in resultados]
    at_contents = [r['at_content'] for r in resultados]
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 15))
    
    ax1.bar(nomes, tamanhos)
    ax1.set_title('Tamanho dos cromossomos')
    ax1.set_xlabel('Cromossomos')
    ax1.set_ylabel('Tamanho (pb)')
    ax1.tick_params(axis='x', rotation=45)
    
    ax2.bar(nomes, gc_contents, label='GC')
    ax2.bar(nomes, at_contents, bottom=gc_contents, label='AT')
    ax2.set_title('Composição de bases dos cromossomos')
    ax2.set_xlabel('Cromossomos')
    ax2.set_ylabel('Porcentagem (%)')
    ax2.tick_params(axis='x', rotation=45)
    ax2.legend()
    
    ax3.scatter(gc_contents, tamanhos)
    for i, nome in enumerate(nomes):
        ax3.annotate(nome, (gc_contents[i], tamanhos[i]))
    ax3.set_title('Relação entre conteúdo GC e tamanho dos cromossomos')
    ax3.set_xlabel('Conteúdo GC (%)')
    ax3.set_ylabel('Tamanho (pb)')
    
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
        
        print("\nResultados detalhados:")
        for r in resultados:
            print(f"{r['nome']}:")
            print(f"  Tamanho: {r['tamanho']} pb")
            print(f"  Conteúdo GC: {r['gc_content']:.2f}%")
            print(f"  Conteúdo AT: {r['at_content']:.2f}%")
            print()
