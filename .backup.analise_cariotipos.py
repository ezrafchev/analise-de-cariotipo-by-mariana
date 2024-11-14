
import os
import re
from Bio import SeqIO
from Bio.SeqUtils import GC, molecular_weight
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

def analisar_cariotipos(arquivo_fasta):
    sequencias = list(SeqIO.parse(arquivo_fasta, "fasta"))
    resultados = []
    
    for seq in sequencias:
        tamanho = len(seq)
        gc_content = GC(seq.seq)
        at_content = 100 - gc_content
        peso_molecular = molecular_weight(seq.seq)
        regioes_repetitivas = identificar_regioes_repetitivas(seq.seq)
        densidade_genes = simular_densidade_genes(tamanho)
        ilhas_cpg = identificar_ilhas_cpg(seq.seq)
        complexidade = calcular_complexidade_sequencia(seq.seq)
        
        resultados.append({
            'nome': seq.id,
            'tamanho': tamanho,
            'gc_content': gc_content,
            'at_content': at_content,
            'peso_molecular': peso_molecular,
            'regioes_repetitivas': regioes_repetitivas,
            'densidade_genes': densidade_genes,
            'ilhas_cpg': ilhas_cpg,
            'complexidade': complexidade
        })
    
    return resultados

def identificar_regioes_repetitivas(seq, min_repeat_length=10):
    repeats = []
    for i in range(len(seq) - min_repeat_length + 1):
        for j in range(i + min_repeat_length, len(seq) + 1):
            subseq = seq[i:j]
            if seq.count(subseq) > 1:
                repeats.append((i, j, len(subseq), seq.count(subseq)))
    return repeats

def simular_densidade_genes(tamanho, gene_density=0.01):
    num_genes = int(tamanho * gene_density)
    gene_positions = sorted(np.random.choice(tamanho, num_genes, replace=False))
    return gene_positions

def identificar_ilhas_cpg(seq, min_length=200, min_gc=50, min_obs_exp=0.6):
    ilhas = []
    for i in range(len(seq) - min_length + 1):
        subseq = seq[i:i+min_length]
        gc_content = GC(subseq)
        obs_exp = calcular_obs_exp_cpg(subseq)
        if gc_content >= min_gc and obs_exp >= min_obs_exp:
            ilhas.append((i, i+min_length))
    return ilhas

def calcular_obs_exp_cpg(seq):
    c_count = seq.count('C')
    g_count = seq.count('G')
    cg_count = seq.count('CG')
    if c_count > 0 and g_count > 0:
        exp_cg = (c_count * g_count) / len(seq)
        return cg_count / exp_cg
    return 0

def calcular_complexidade_sequencia(seq, window_size=100):
    complexidades = []
    for i in range(0, len(seq) - window_size + 1, window_size):
        subseq = seq[i:i+window_size]
        complexidade = len(set(subseq)) / window_size
        complexidades.append(complexidade)
    return np.mean(complexidades)

def visualizar_resultados(resultados):
    nomes = [r['nome'] for r in resultados]
    tamanhos = [r['tamanho'] for r in resultados]
    gc_contents = [r['gc_content'] for r in resultados]
    at_contents = [r['at_content'] for r in resultados]
    
    fig, axes = plt.subplots(3, 2, figsize=(20, 30))
    
    # Tamanho dos cromossomos
    axes[0, 0].bar(nomes, tamanhos)
    axes[0, 0].set_title('Tamanho dos cromossomos')
    axes[0, 0].set_xlabel('Cromossomos')
    axes[0, 0].set_ylabel('Tamanho (pb)')
    axes[0, 0].tick_params(axis='x', rotation=45)
    
    # Composição de bases
    axes[0, 1].bar(nomes, gc_contents, label='GC')
    axes[0, 1].bar(nomes, at_contents, bottom=gc_contents, label='AT')
    axes[0, 1].set_title('Composição de bases dos cromossomos')
    axes[0, 1].set_xlabel('Cromossomos')
    axes[0, 1].set_ylabel('Porcentagem (%)')
    axes[0, 1].tick_params(axis='x', rotation=45)
    axes[0, 1].legend()
    
    # Relação entre conteúdo GC e tamanho
    axes[1, 0].scatter(gc_contents, tamanhos)
    for i, nome in enumerate(nomes):
        axes[1, 0].annotate(nome, (gc_contents[i], tamanhos[i]))
    axes[1, 0].set_title('Relação entre conteúdo GC e tamanho dos cromossomos')
    axes[1, 0].set_xlabel('Conteúdo GC (%)')
    axes[1, 0].set_ylabel('Tamanho (pb)')
    
    # Densidade de genes
    gene_densities = [len(r['densidade_genes']) / r['tamanho'] for r in resultados]
    axes[1, 1].bar(nomes, gene_densities)
    axes[1, 1].set_title('Densidade de genes (simulada)')
    axes[1, 1].set_xlabel('Cromossomos')
    axes[1, 1].set_ylabel('Genes por base')
    axes[1, 1].tick_params(axis='x', rotation=45)
    
    # Complexidade da sequência
    complexidades = [r['complexidade'] for r in resultados]
    axes[2, 0].bar(nomes, complexidades)
    axes[2, 0].set_title('Complexidade da sequência')
    axes[2, 0].set_xlabel('Cromossomos')
    axes[2, 0].set_ylabel('Complexidade média')
    axes[2, 0].tick_params(axis='x', rotation=45)
    
    # Ideograma do cariótipo
    ax_ideogram = axes[2, 1]
    y_positions = np.arange(len(nomes))
    ax_ideogram.barh(y_positions, tamanhos, height=0.5)
    ax_ideogram.set_yticks(y_positions)
    ax_ideogram.set_yticklabels(nomes)
    ax_ideogram.set_title('Ideograma do Cariótipo')
    ax_ideogram.set_xlabel('Tamanho (pb)')
    
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
            print(f"  Peso molecular: {r['peso_molecular']:.2f}")
            print(f"  Número de regiões repetitivas: {len(r['regioes_repetitivas'])}")
            print(f"  Número de genes (simulado): {len(r['densidade_genes'])}")
            print(f"  Número de ilhas CpG: {len(r['ilhas_cpg'])}")
            print(f"  Complexidade da sequência: {r['complexidade']:.4f}")
            print()
