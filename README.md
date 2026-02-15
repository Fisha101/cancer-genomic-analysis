# 🧬 Cancer Genomic Analysis  

This repository contains a Python-based computational genomics project focused on sequence-level characterization of selected human oncogenes.  

The project demonstrates foundational bioinformatics concepts including sequence parsing, compositional analysis, and open reading frame (ORF) detection implemented without external bioinformatics libraries.

---

## 📘 Project Overview  

This analysis examines genomic features of key cancer-related genes by:

- Parsing DNA sequences from FASTA files  
- Calculating GC content (%)  
- Scanning all three forward reading frames  
- Identifying the longest Open Reading Frame (ATG → in-frame stop codon)  
- Generating automated CSV reports  
- Producing GC content visualization  

The goal is to demonstrate algorithmic understanding of coding structure and sequence composition in human oncogenes.

---

## 🧬 Genes Analyzed  

- **TP53** – tumor suppressor gene involved in DNA damage response  
- **KRAS** – proto-oncogene regulating cell proliferation signaling  
- **MYC** – transcription factor driving cell cycle progression  
- **BRCA1** – DNA repair–associated tumor suppressor gene  

---

## 🧠 Repository Files  

### 🔹 Code  
- [`main.py`](main.py) – Complete analysis pipeline: FASTA parsing, GC-content calculation, multi-frame ORF detection, CSV export, and visualization.

### 🔹 Input Data  
- [`data/`](data/) – Directory containing FASTA files of selected human oncogenes.

### 🔹 Output  
- [`report.csv`](report.csv) – Generated genomic metrics for each analyzed gene.  
- [`plots/gc_content.png`](plots/gc_content.png) – Bar chart of GC content across genes.

---

## 📈 Interpretation  

Preliminary results show variation in GC content across oncogenes, reflecting gene-specific sequence composition patterns.  

The identified longest ORFs correspond to substantial coding regions consistent with functional protein-coding genes, confirming correct reading-frame logic and biologically meaningful ORF detection.

---

## ▶️ How to Run  

    python3 main.py
