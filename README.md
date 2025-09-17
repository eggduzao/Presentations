# 🌟 Presentations: Codes, Data, and Notes 🌟  

Welcome to **Presentations**, a repository designed as a living archive of all my **presentation codes, scripts, datasets, visualizations, and lecture notes**.  
This repo is not just a code dump — it’s a curated, well-documented space to keep track of everything I’ve shown, taught, or explained in talks, tutorials, workshops, and academic presentations.  

Whether you are here to **reuse snippets**, **explore genomics demonstrations**, or simply browse the structure of past presentations, this repository is intended to be:  

- 📚 **Educational**: All content is written with clarity, with pointers for self-learning.  
- 🔬 **Reproducible**: Code and data references are stored as much as possible to allow direct re-use.  
- 🎨 **Polished**: Each presentation is not only scientific but also formatted for readability and elegance.  
- 🔗 **Connected**: All entries are hyperlinked with internal anchors for smooth navigation.  

---

## ✨ Features of this Repository  

✔️ **Organized by Presentation Date** (chronological order).  
✔️ **Code + Notes + Data pointers** included in each folder.  
✔️ **Table of Contents (TOC)** maintained here for quick reference.  
✔️ **Direct Anchors** → Jump to detailed sections below.  
✔️ **MIT Licensed** → Free to use, modify, adapt, and remix.  

---

## 📑 Table of Contents  

Each row corresponds to a specific presentation. The table will expand as new entries are added, sorted by **increasing presentation day**.  

| Date       | Presentation Name | Description | Link |
|------------|------------------|-------------|------|
| 2025-09-17 | IGV Tutorial     | A guided demonstration of genome visualization in IGV, exploring peaks, tracks, and common pitfalls. | [Details](#2025-09-17-igv-tutorial) |

---

## 🗂️ Presentation Details

### 2025-09-17: IGV Tutorial  
📂 Path: `<Presentations>/ReGENE/2025_09_17_Artigo_IGV/`  

This presentation documents my exploration of **IGV (Integrative Genomics Viewer)** — not only as a visualization tool, but as a narrative of how subtle technical details, file formats, and genome builds can radically alter what we see (or fail to see).  
It is meant to serve both as a **technical walkthrough** and a **cautionary tale** for anyone preparing IGV demos: things that look trivial (colors, alignments, SNPs) may actually be controlled by hidden metadata, genome builds, or file preprocessing.  

#### Key Themes and Lessons  

- **BAM File Behavior in IGV**
  - Default behavior: BAMs load with colors only for mismatches.
  - Exception case: 1000 Genomes BAMs sometimes load fully colored.  
    → Root cause: **SAM header read-group differences** and **genome build mismatches** (hg19 vs hg38).  
  - Lesson: Always check `@RG` headers and confirm genome build before assuming IGV is “misbehaving.”

- **VCF Handling and Indexing**
  - Attempted to index `ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` with `tabix`.  
    - Standard `tabix -p vcf` failed because the file is bgzipped differently.  
    - Workaround: `gunzip -c | bgzip -c` rewrapping.  
  - Key point: **1000 Genomes "ALL" VCFs often require rewrapping before indexing**.

- **Genome Build Mismatches**
  - Initially loaded hg38 genome in IGV with 1000 Genomes BAMs (which are aligned to hg19).  
  - Result: *every single base colored*, looking like a systematic variant at every position.  
  - Core insight: **if every base is highlighted, suspect a reference mismatch** rather than biology.  

- **“What’s Interesting to Show?”**
  - Realization: many regions look identical or unremarkable when naively browsing BAMs.  
  - Suggested strategies:  
    - Use known variant loci (from dbSNP / 1000 Genomes VCFs).  
    - Demonstrate promoter/enhancer divergence in K562 H3K4me3 data.  
    - Show tissue-specific alternative splicing (e.g., SLC25A3 with heart vs liver RNA-seq).  

- **ENCODE Dataset Navigation**
  - Differences between “alignments” and “transcriptome alignments”:  
    - *Alignments*: reads mapped to the full reference genome.  
    - *Transcriptome alignments*: reads mapped only to transcript models, useful for isoform expression quantification.  
  - Duplicate-looking BAMs of the same replicate with different sizes: often represent different processing versions (filtered vs unfiltered, different aligners, or inclusion/exclusion of multimapping reads).

- **Comparisons Across Protocols**
  - PolyA RNA-seq vs RAMPAGE:  
    - **Not directly comparable** — PolyA captures bulk mRNA, RAMPAGE is designed for transcription start sites (TSS).  
    - Both can be shown side by side, but with the caveat that they answer different biological questions.

- **Tissue vs Cell Type Metadata**
  - ENCODE “heart tissue” or “liver tissue” means *bulk tissue RNA-seq*, not sorted cell types.  
  - Important caveat: results represent **mixed populations of cells** rather than purified cell-type signals.  

#### Narrative Arc for the Presentation  

1. **Setting the Scene**: Load a BAM, expect mismatches only. Suddenly — everything is colored. Suspicion arises.  
2. **Diagnosis Journey**: Realize that hidden SAM headers and reference genome mismatches explain the chaos.  
3. **Broader Exploration**: Attempt to connect BAMs with VCFs, rewrap files, and make sense of what’s worth showing.  
4. **The “Showpiece” Examples**:  
   - K562 H3K4me3 promoter divergence.  
   - SLC25A3 alternative splicing in heart vs liver RNA-seq.  
5. **Take-Home Message**:  
   - IGV is powerful, but **only as good as your metadata discipline**.  
   - Always document genome builds, file preparation steps, and dataset provenance before presenting.  

🚀 This section of the repository thus represents both the **technical playbook** and the **conceptual reflections** that emerged while wrestling with IGV.  

---

---

## 📜 License  

This repository is licensed under the **MIT License**.  
You are free to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the material, provided that proper attribution is maintained.  

---

## 🤝 Contributing  

Contributions are welcome! If you’ve spotted a typo, have a suggestion for clearer explanations, or want to add complementary code or datasets, feel free to open a Pull Request.  

---

## 🌐 Stay Connected  

- ✉️ For questions or discussions, open an Issue in this repository.  
- 🔬 More about me: *[Add your personal website / academic page here]*  
- 🧬 Future topics: Expect more on genomics, bioinformatics pipelines, machine learning, and teaching demos.  

---

💡 *This repository is meant to grow — just like knowledge itself. Each presentation tells a story, and together they form a timeline of learning and sharing.*  