import os
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_fasta(path):
    sequences = {}
    current_gene = None
    current_sequence = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_gene is not None:
                    sequences[current_gene] = "".join(current_sequence)

                current_gene = line[1:].strip()
                current_sequence = []
            else:
                current_sequence.append(line)

    if current_gene is not None:
        sequences[current_gene] = "".join(current_sequence)

    return sequences


def kmer_frequencies(seq, k):
    counts = Counter()

    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        counts[kmer] += 1

    total = sum(counts.values())

    if total == 0:
        return {}

    freqs = {}
    for kmer, count in counts.items():
        freqs[kmer] = count / total

    return freqs

def cosine_similarity_matrix(X):
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    Xn = X / norms
    return Xn @ Xn.T


def run_kmer_pipeline():
    data_folder = "data"
    all_sequences = {}

    for filename in os.listdir(data_folder):
        if filename.endswith(".fasta"):
            path = os.path.join(data_folder, filename)
            sequences = parse_fasta(path)
            all_sequences.update(sequences)

    print("Total genes loaded:", len(all_sequences))
    print("Gene names:", list(all_sequences.keys()))

    # ---- build gene_to_freq (k-mer fingerprints) ----
    k = 4  # change if needed
    gene_to_freq = {}

    for gene, seq in all_sequences.items():
        seq = seq.upper().replace(" ", "").replace("\n", "")
        gene_to_freq[gene] = kmer_frequencies(seq, k)

    print("Genes processed:", len(gene_to_freq))

    # ---- build feature matrix ----
    out_folder = "results/kmer_fingerprints"
    os.makedirs(out_folder, exist_ok=True)
    os.makedirs(os.path.join(out_folder, "figures"), exist_ok=True)

    all_kmers = set()
    for g, freq in gene_to_freq.items():
        all_kmers.update(freq.keys())
    all_kmers = sorted(all_kmers)

    genes = list(gene_to_freq.keys())
    matrix = []
    for g in genes:
        row = [gene_to_freq[g].get(km, 0.0) for km in all_kmers]
        matrix.append(row)

    df_features = pd.DataFrame(matrix, index=genes, columns=all_kmers)
    df_features.to_csv(os.path.join(out_folder, f"kmer_features_k{k}.csv"))

    # ---- cosine similarity ----
    X = df_features.values.astype(float)
    sim = cosine_similarity_matrix(X)
    df_sim = pd.DataFrame(sim, index=genes, columns=genes)
    df_sim.to_csv(os.path.join(out_folder, f"kmer_cosine_similarity_k{k}.csv"))

    # ---- top-3 neighbors ----
    rows = []
    for g in genes:
        s = df_sim.loc[g].drop(index=g).sort_values(ascending=False)
        top = s.head(3)
        rows.append({
            "gene": g,
            "neighbor_1": top.index[0] if len(top) > 0 else "",
            "sim_1": float(top.iloc[0]) if len(top) > 0 else np.nan,
            "neighbor_2": top.index[1] if len(top) > 1 else "",
            "sim_2": float(top.iloc[1]) if len(top) > 1 else np.nan,
            "neighbor_3": top.index[2] if len(top) > 2 else "",
            "sim_3": float(top.iloc[2]) if len(top) > 2 else np.nan,
        })

    df_top = pd.DataFrame(rows)
    df_top.to_csv(os.path.join(out_folder, f"top_neighbors_k{k}.csv"), index=False)

    # ---- heatmap ----
    plt.figure(figsize=(8, 6))
    plt.imshow(df_sim.values, aspect="auto")
    plt.colorbar(label="Cosine similarity")
    plt.xticks(range(len(genes)), genes, rotation=90, fontsize=7)
    plt.yticks(range(len(genes)), genes, fontsize=7)
    plt.title(f"k-mer similarity heatmap (k={k})")
    plt.tight_layout()
    plt.savefig(os.path.join(out_folder, "figures", f"kmer_similarity_heatmap_k{k}.png"), dpi=200)
    plt.close()

    print("Saved results to:", out_folder)

    return gene_to_freq


gene_to_freq = run_kmer_pipeline()






   
    
