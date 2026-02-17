# 🧬 Cancer Genomic Analysis

This repository contains two complementary computational genomics pipelines focused on sequence-level characterization and similarity analysis of selected human oncogenes.

Both modules are implemented in Python using core scientific libraries (NumPy, Pandas, Matplotlib) without external bioinformatics frameworks.

---

# 📘 Project 1 — Genomic Structure Analysis

This module focuses on sequence composition and coding structure.

## 🔬 Methods

- FASTA parsing  
- GC content calculation (%)  
- Three forward reading frame scanning  
- Longest Open Reading Frame (ATG → in-frame stop codon) detection  
- Automated CSV report generation  
- GC content visualization  

## ▶️ Run

```bash
python3 main.py
```

## 📤 Output

- `report.csv` — GC content and longest ORF statistics per gene  
- `plots/gc_content.png` — GC content comparison plot  

## 🧠 Biological Insight

Variation in GC content reflects gene-specific compositional patterns.

Detected longest ORFs correspond to substantial coding regions consistent with functional protein-coding genes, validating reading-frame logic and ORF detection implementation.

---

# 📘 Project 2 — k-mer Genomic Fingerprints & Similarity

This module represents each gene as a normalized k-mer frequency vector and performs compositional similarity analysis.

## 🔬 Methods

- k-mer frequency computation (default k = 4)  
- Feature matrix construction (genes × kmers)  
- Cosine similarity calculation  
- Top-3 nearest neighbor identification  
- Similarity heatmap visualization  

## ▶️ Run

```bash
python3 modules/kmer_fingerprints/kmer.py
```

## 📤 Output

All results are saved to:

```
results/kmer_fingerprints/
```

Generated files:

- `kmer_features_k4.csv` — feature matrix  
- `kmer_cosine_similarity_k4.csv` — similarity matrix  
- `top_neighbors_k4.csv` — top-3 most similar genes  
- `figures/kmer_similarity_heatmap_k4.png` — similarity heatmap  

## 🧠 Computational Insight

Cosine similarity compares directional patterns of normalized k-mer frequency vectors, enabling compositional similarity analysis independent of gene length.

This approach transforms nucleotide sequences into quantitative feature representations suitable for downstream computational analysis.

---

# 🧬 Genes Analyzed

- **TP53** — tumor suppressor involved in DNA damage response  
- **KRAS** — proto-oncogene regulating proliferation signaling  
- **MYC** — transcription factor driving cell cycle progression  
- **BRCA1** — DNA repair–associated tumor suppressor  

---

# 🗂 Repository Structure

```
.
├── main.py
├── modules/
│   └── kmer_fingerprints/
│       └── kmer.py
├── data/
├── plots/
├── results/
│   └── kmer_fingerprints/
└── README.md
```

---

# ⚙️ Requirements

- Python 3.9+
- numpy
- pandas
- matplotlib

Install dependencies if needed:

```bash
pip install numpy pandas matplotlib
```

---

# 🎯 Project Objective

To demonstrate foundational computational genomics skills including:

- sequence parsing  
- compositional analysis  
- reading-frame logic  
- feature engineering  
- vector-based similarity analysis  
- reproducible reporting  
- scientific visualization  

The combined pipelines provide a structured framework for sequence-level characterization of oncogenes using pure Python.

