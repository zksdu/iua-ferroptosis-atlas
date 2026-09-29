# Cover Letter — Paper 2（2026-09-29 v1.1 robustness edition）

**To**: Editor-in-Chief, Frontiers in Immunology
**Manuscript**: Ferroptosis defense collapse, not activation, defines epithelial susceptibility in intrauterine adhesions: a dual-resolution transcriptomic atlas
**Article type**: Original Research

Dear Editor,

Intrauterine adhesions (IUA, Asherman syndrome) are a leading yet poorly understood cause of infertility. A recent landmark study (Zhu et al., Free Radic Biol Med 2023) established that ferroptosis participates in endometrial fibrosis and proposed ferroptosis induction in epithelium as the pathogenic event. We asked a simple but decisive question at single-cell resolution: **is the IUA endometrium executing ferroptosis — or has it lost the defenses that prevent it?**

Using bulk RNA-seq of severe IUA (GSE224093, 7 v 7), a single-cell atlas of the Asherman endometrium (GSE215968, 106,400 cells, 17 cell types), an independent within-patient direction cohort (GSE160633), and paired pre/post-therapy biopsies (GSE311899), we find that:

1. Bulk ferroptosis pathway scores are **lower** — not higher — in IUA, with concordant collapse of the GPX4/DHODH/KEAP1/NCOA4 defense axis and Wnt reprogramming (ferroptosis–Wnt ρ=0.859);
2. Epithelium is the **only major cell type with concordant dual deterioration** (activity +0.145, defense −0.244), most severe in the disease-specific AS-epithelium subtype (−0.335) — reproducing, at population scale, the epithelial focus of Zhu et al. while reinterpreting it as a **susceptible state**;
3. **Macrophages are the iron sink** of the fibrotic endometrium (only significantly elevated iron load; highest absolute level), with ferroptosis–inflammation coupling **myeloid-anchored** (DC 0.31 > macrophages 0.30 >> epithelium 0.07) — a niche signature we observe in parallel analyses of preeclamptic placenta, suggesting a shared myeloid-anchored iron niche across fibrotic pregnancy disease;
4. Bulk-level HMGB1 decrease is fully explained by **compartmentalized redistribution**: parenchymal collapse concurrent with myeloid activation (DC +0.291);
5. Findings survive sample-level re-analysis (7/7 direction-consistent; epithelial defense collapse p=0.0067, Hedges g=−1.21);
6. a dedicated robustness suite — leave-one-out gene-set perturbation (52/52 runs direction-consistent across the four module panels), independent KEGG ferroptosis gene-set cross-scoring, per-cell-type Wnt coupling localization, myeloid M1/M2 subtyping, and bulk composition-adjusted re-testing — addresses gene-set bias, cell-composition confounding and alternative hypotheses without any new experiments (Supplementary Figures S3–S6, Table S6).

**Why this matters**: the distinction is therapeutic. If ferroptosis were already executed, inhibitors would be too late; a collapsed defense state that is **potentially reversible** — reversibility currently supported by murine ferrostatin-1 evidence and trial-reported endometrial regeneration, while paired post-therapy biopsies (GSE311899) do not yet resolve epithelial defense recovery at mRNA level, which we report transparently — argues for **restoring epithelial anti-ferroptotic defenses** (e.g., GPX4/System xc− axis), and nominates the myeloid iron niche as a targetable microenvironmental component.

We believe this work fits Frontiers in Immunology's readership: it connects ferroptosis, innate myeloid biology and tissue fibrosis, and provides a cell-type-resolved framework that reframes an active therapeutic debate.

All analyses are fully reproducible: the complete analysis code, intermediate result tables, and publication figures have been assembled for release in a public code repository upon submission.

**Declarations**: all data are from public GEO repositories; original analyses; no conflicts; all authors have approved submission. The manuscript is not under consideration elsewhere.

Sincerely,
Bing Song
The Third Affiliated Hospital of Guangzhou Medical University
bingsong2012683034@gzhmu.edu.cn

---

## 投稿要点备忘（内部，不随信发出）
- 审稿人 1 预判：GSE215968 对照为 WOI（月经周期混杂）→ 已在 Limitations 声明，可补周期分层敏感性（若被要求）
- 审稿人 2 预判：mRNA 级证据 → 回应：FRBM 2023 已提供蛋白级 IHC（4-HNE/GPX4），我们与其蛋白观察兼容；补充声明未来湿实验验证
- 审稿人 3 预判：GSE311899 同组来源 + 无健康对照 → 已在 3.7 与 Methods 明确 "additional paired cohort from the same trial"，不称独立
- 与 Zhu et al. 关系定位：尊重+延伸（"their epithelial observation is correct and now explained"），避免对抗性措辞
