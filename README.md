# 🧬 Sjögren-Larsson Syndrome (SLS) VUS Analysis

This project investigates **variants of unknown significance (VUS)** in the **ALDH3A2 gene**, which is associated with **Sjögren-Larsson Syndrome (SLS)**—a rare autosomal recessive neurocutaneous disorder. The analysis combines **evolutionary biology**, **phylogenetics**, and **bioinformatics tools** to assess the likelihood of these VUS being pathogenic or benign.

---

## 📚 Introduction

Sjögren-Larsson Syndrome (SLS) is a rare genetic condition first described in Sweden in 1957. It affects **muscle tone**, **skin**, **retinal pigmentation**, and **intellectual development**, caused by **mutations in the ALDH3A2 gene** on chromosome 17. This gene encodes the fatty aldehyde dehydrogenase (FALDH) enzyme.

This project focuses on the **evolutionary history of ALDH3A2** and characterization of **variants of uncertain significance (VUS)** using **phylogenetic trees**, **sequence alignment**, and **SIFT scores**.

---

## 🔬 Objective

- Analyze ALDH3A2 orthologs and paralogs
- Visualize evolutionary relationships with **Maximum Likelihood trees**
- Assess the **tolerability** of VUS based on amino acid conservation and phylogenetic distance
- Interpret functional consequences using **SIFT** and amino acid frequency analysis

---

## 🧪 Tools & Methodology

- **Databases**:
  - [Ensembl](https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000072210)
  - [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/)
  - [LOVD ALDH3A2](https://databases.lovd.nl/shared/genes/ALDH3A2)
  - NORD Rare Diseases

- **Software & Libraries**:
  - MEGA11 (Muscle, Maximum Likelihood Trees, Bootstrapping)
  - Python (sequence mutation recreation, FASTA processing)
  - SIFT (tolerability predictions)
  - FigTree (rerooting and visualizing NJ trees)
  - Pylab (visual analysis)

- **Approach**:
  - Download orthologous and paralogous sequences of ALDH3A2
  - Perform multiple sequence alignments
  - Construct phylogenetic trees using MEGA11
  - Analyze amino acid frequencies and tolerance scores
  - Evaluate VUS using tree positioning and SIFT scores

---

## 📈 Results Summary

- **Phylogenetic Trees** (Figures 1–6): Constructed using Maximum Likelihood and Neighbor Joining methods; annotated with bootstrap values and visual color-coding for variant presence.

![image](https://github.com/user-attachments/assets/6d14d698-9d4c-4acf-abd3-2d38f20b516c)

![image](https://github.com/user-attachments/assets/b38707eb-bcc3-48d4-bc42-32a5466aa41a)

![image](https://github.com/user-attachments/assets/9952c3a2-bb55-42ef-853c-93fd1dd68403)

- 

- **VUS Interpretation**:
  - Variants close to mammalian lineages = more tolerable
  - Distant sequences (e.g., yellow color in rerooted tree) suggest intolerant or deleterious variants
 
    ![image](https://github.com/user-attachments/assets/fdb16a5a-9bb7-4ba9-9864-b2c5bd3c9632)

- **Tables**:
  - **Table 1**: Amino acid frequency per position
  - **Table 2**: Tolerance scores of each VUS

> Example Insight: A VUS at position 339 (I → V) may be tolerable; while the same change at position 352 is less clear and requires further validation.

---

## 🧾 Files Included

- `aligned_sequences.fasta` — Aligned orthologs/paralogs
- `vus_variants.csv` — Annotated VUS with positions and substitutions
- `mega_trees/` — Exported `.nwk` tree files from MEGA
- `figures/` — Visualized trees from MEGA and FigTree
- `tolerance_scores.csv` — Scores for each VUS
- Python scripts:
  - `recreate_variants.py`
  - `blast_merge_clean.py`
  - `analyze_frequencies.py`

---

## 📌 Conclusion

This project demonstrates how **phylogenetic analysis** and **bioinformatic predictions** can help classify **VUS** in rare diseases such as SLS. Using tools like **MEGA**, **SIFT**, and evolutionary conservation, it is possible to infer functional relevance of VUS and support variant interpretation efforts in clinical genomics.

---


## 📂 License

This project is for academic and research purposes only.

