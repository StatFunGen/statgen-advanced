# statgen-advanced

These notes are for trainees with quantitative backgrounds but without formal training in statistical genetics, who have encountered these methods in the literature but have not yet worked with them hands-on. For statisticians wanting to catch up on genetics applications, these notes also provide essential conceptual foundations and key assumptions that geneticists make when modeling data.

These notes are not organized by method, by paper, or by software tool. Instead, we organize by scientific question. For each question, we focus on what problem we are trying to solve, what assumptions we are making, and what generative model most naturally describes how the data arise. Once these foundations are clear, existing methods become natural solutions, and their limitations become obvious.

Think of it like building a Lego model to represent something in the real world. The statistical building blocks (likelihoods, priors, latent variables, hierarchical structures) are the pieces available. Our goal is to focus on designing the blueprint that captures the essential features of the biological reality, while keeping the available blocks in mind. The details of assembling specific kits will inevitably be discussed, but they are not the focus. When one understands what the design requires and what connections matter, one will know how to select and combine blocks to satisfy those requirements. With this foundation, one can read new methods papers and recognize the same underlying ideas, and feel comfortable adapting or extending existing approaches for new problems (and only worry about various tricky details at this point).

As an example, consider allele-specific expression (ASE) QTL analysis. Total expression reflects the sum of transcripts from both haplotypes; ASE measures their difference within heterozygotes. The same genetic effect parameter underlies both, appearing as dosage effect $(0, 1, 2)$ in total expression and haplotype difference $(-1, 0, +1)$ in ASE. Because sum and difference are conditionally independent, ASE adds information about genetic effects beyond total expression from the same samples, effectively increasing sample size. The within-individual comparison also cancels individual-level confounders (which affect both haplotypes equally), and haplotype difference in ASE provides different correlation (LD) patterns than conventional genotype dosage, thus improving fine-mapping resolution. These advantages motivate incorporating ASE into QTL analysis. [RASQUAL](https://www.nature.com/articles/ng.3467) implemented a rigorous generative model with Negative Binomial total counts and Beta-Binomial allele-specific counts sharing genetic effect parameters; [mixQTL](https://www.nature.com/articles/s41467-021-21592-8) later achieved scalability through Gaussian approximations and [WASP](https://github.com/bmvdgeijn/WASP) preprocessing, trading some modeling rigor for computational efficiency suitable for large-scale analysis. One can extend this framework further by adding local ancestry modeling and fine-mapping under different likelihoods, following the same approach to motivation and generative modeling.

## Overview of Topics

These notes organize into five themes. 

The first three themes represent our understanding of how genetic association relates to biological questions. When studying genetic associations, we fundamentally care about two things: mapping specific effects to variants or contexts (Theme 1) or predicting aggregate outcomes from polygenic architecture (Theme 2). Theme 3 takes both perspectives further by adding causal inference: we can ask mapping questions (which gene's expression causes disease?) or prediction questions (can we predict disease risk through genetically predicted expression?). Causality is necessarily narrower than association—it requires stronger assumptions—but builds on the same conceptual foundations. 

The last two address how we adapt our models to specific data types or practical computational constraints. Theme 4 examines how biological data-generating processes inform model design. Unlike simple GWAS phenotypes, molecular data (RNA-seq counts, methylation proportions, splicing ratios) arise from well-understood biological mechanisms that suggest specific distributional choices and model structures. Theme 5 addresses strategic simplification of rigorous generative models for computational tractability. 

Throughout, the same building blocks, mostly introduced in our ["statgen-primer" notes](https://statfungen.github.io/statgen-primer), appear in different combinations depending on the scientific question.

### Theme 1: Mapping of shared vs. specific effects

The core question in genetic mapping is: which variants are associated with the trait, and how do effects vary across contexts? We ask these questions whether comparing the same trait across studies (meta-analysis), the same trait across ancestries (cross-ancestry analysis), or different traits at the same locus (colocalization).
 
In sum, both genetic effects and confounders can be shared or context-specific. In meta-analysis, we typically assume effects are shared (or allow for heterogeneity) while residuals are study-specific. In cross-ancestry fine-mapping, effects may be shared but LD patterns (a type of confounder) differ by ancestry. In multi-tissue QTL analysis, confounders may be shared across tissues belonging to the same biological sample donor, while effects vary by cell or tissue types even within the same donor.

Fine-mapping is discussed here because LD is simply a particular type of confounding, where correlated variants make it hard to identify causal ones. The challenge is high-dimensional (many correlated variants) and requires variable selection because we don't know which effects are real. Annotation priors for fine-mapping are generative model details that help distinguish confounders from true signals, not a separate concept.

Note that colocalization asks whether two traits share a causal variant at a locus. This is a mapping question (Theme 1), not a causal inference question (Theme 3). Colocalization results are biologically suggestive of shared mechanism, but the statistical model makes no causal assumptions about one trait affecting another.

| Scientific Question | Core Concepts | Methods/Tools |
|---------------------|---------------|---------------|
| Which variants are likely causal given LD? | Fine-mapping, credible sets, variable selection, functional priors | SuSiE, FINEMAP, PolyFun |
| How do we combine evidence across studies? | Effect sharing and heterogeneity, meta-analysis as special case of multivariate modeling | METAL, METASOFT, mtag, mashr|
| Do two traits share causal variants? | Colocalization, shared genetic architecture | coloc, SuSiE-coloc, ColocBoost |
| How do we leverage LD diversity across populations? | Cross-ancestry fine-mapping, ancestry-specific LD | SuSiEx, MESuSiE, MultiSuSiE, SuShiE |
| How do genetic effects vary across contexts? | Multi-trait GWAS, multi-context QTL, effect heterogeneity | mtCOJO, mvSuSiE |

### Theme 2: Prediction with polygenic architecture

Sometimes we don't care which specific variants are causal. Instead, we want to understand aggregate properties: How much of trait variation is genetic? Which tissues or cell types are enriched? Can we predict phenotype from genotype? This is the polygenic perspective, where the signal that cannot be mapped to individual variants still contributes meaningfully to prediction and biological insight.

In contrast to theme 1 where fine-mapping asks "which?", here heritability asks "how much?" PRS asks "what can we predict?". LDSC leverages LD structure not to identify causal variants but to estimate total heritability and its enrichment across functional categories. The unmappable portion of genetic signal, which would be noise in Theme 1, becomes useful information in Theme 2.

MAGMA fits here as gene-level aggregation: it tests whether variants near a gene collectively associate with a trait, without invoking causal assumptions about gene expression affecting the trait (that would be Theme 3). Genetic correlation between traits also typically belongs here, as we usually refer to genome-wide correlation computed in a polygenic framework.

| Scientific Question | Core Concepts | Methods/Tools |
|---------------------|---------------|---------------|
| How much trait variation is genetic? | Heritability, variance components | GREML, GCTA, LDSC, HESS, LDAK, HDL |
| Are genetic effects correlated across traits? | Genetic correlation, pleiotropy | bi-LDSC, Popcorn, HESS |
| Which tissues or annotations are enriched? | Heritability partitioning, functional enrichment | S-LDSC |
| Can we predict phenotype from genotype? | Polygenic risk scores, prediction accuracy | PRSice, LDpred, PRS-CS, PRS-CSx |

### Theme 3: Causal chain from molecular phenotype to disease

This theme applies the concepts from Themes 1 and 2 most typically to the relationship between molecular phenotypes and disease, with additional statistical assumptions that enable causal inference. The key distinction is the instrumental variable framework: we use genetic variants as instruments to test whether an exposure (e.g., gene expression) causally affects an outcome (e.g., disease risk).

It is important to clarify that "causal" here refers to the statistical modeling framework, not biological mechanism. MR and TWAS both rely on instrumental variable assumptions (relevance, independence, exclusion restriction), and all TWAS methods can be viewed as two-sample MR with gene expression as the exposure. The differences are practical (expression prediction, correlated instruments) rather than conceptual. This unification helps tracking what each method assumes and where it might fail, particularly through horizontal pleiotropy. Colocalization, by contrast, asks whether GWAS and QTL share a causal variant but makes no statistical claim about expression causing disease.

| Scientific Question | Core Concepts | Methods/Tools |
|---------------------|---------------|---------------|
| Does gene expression causally affect disease? | TWAS as MR, instrumental variables, predicted expression | PrediXcan, FUSION, MultiXcan, mr.mash, CoMM, cTWAS |
| Does exposure X cause outcome Y? | MR assumptions, horizontal pleiotropy, instrument selection | TwoSampleMR, MR-Egger, MR-PRESSO, MRAID |
| Which among multiple correlated exposures is the driver? | Multivariate MR (MVMR), direct vs. indirect effects, pleiotropy correction | MVMR-Egger, MVMR-cML, GRAPPLE, Cis-MRBEE|
| Can we distinguish causality from pleiotropy? | Horizontal pleiotropy testing, robust MR | PMR-Egger, CAUSE |
| How do we integrate QTL and GWAS evidence for causality? | multi-omics MR | SMR, ... |


### Theme 4: Generative models for molecular phenotypes

This theme examine specific generative modeling to different molecular data types. Unlike GWAS where the phenotype is relatively simple (a quantitative trait or case-control status), molecular phenotypes have complex data generating processes that we understand from biology. Building generative models that respect these structures can improve power and interpretation.

This is the Lego analogy in action: we know the biology of RNA-seq count data or bisulfite sequencing, so we can choose distributional assumptions (Negative Binomial, Beta-Binomial) and model structure (haplotype effects, overdispersion) that match reality, followed by tools discussed in previous themes to implement these specific models.

| Scientific Question | Core Concepts | Methods/Tools |
|---------------------|---------------|---------------|
| How do we model splicing variation? | Junction usage, Dirichlet-Multinomial, intron clusters | Leafcutter, sQTLseekeR, ISSAC |
| How do we model methylation? | Beta-distributed outcomes, spatial correlation | smash, fSuSiE |
| How do we model protein abundance? | pQTL mapping, measurement noise, missing data | ... |
| How do we leverage allele-specific information? | ASE, haplotype-aware models, conditional independence | RASQUAL, mixQTL, WASP |
| How do we handle single-cell QTL? | Cell type composition, pseudobulk, mixed effects | ... |

### Theme 5: Scalability and computational approximations



Themes 1-4 emphasize building full generative models that capture biological reality. These models are conceptually complete but often computationally intractable at biobank scale. As genetic datasets grow (hundreds of thousands to millions of individuals), statistical methods often become computationally intractable. Theme 5 addresses how we strategically simplify these generative models through approximations that preserve core conceptual framework while enabling computation. The key is understanding what we're dropping and why it still works.

We emphasize tradeoffs to helps practitioners choose appropriate methods between rigorous generative models and scalable approximations.

Note: Basic linear models, mixed models for relatedness, and population structure correction are covered in [statgen-primer](https://statgen.github.io/statgen-primer/). This course focuses on advanced methods for large-scale data.

| Scientific Question | Core Concepts | Methods/Tools |
|---------------------|---------------|---------------|
| How do we run GWAS on biobank-scale data? | Scalable mixed models, sparse GRM, approximations | BOLT-LMM, SAIGE, REGENIE, fastGWA |
| How do we work with summary statistics? | LD reference panels, avoiding individual-level data | Summary-stat-based methods throughout |

---

**What is NOT discussed**

These notes focuses on statistical genetics methodology. Several related topics are beyond our current scope:

- *Variant effect prediction via AI/ML*: Methods like AlphaMissense, ESM1b, and Enformer use deep learning to predict variant effects from sequence or structure. While increasingly important, these are machine learning predictions without a strong genetics motivation, and are not covered here.
- *Network and pathway analysis*: Methods like DEPICT and pathway enrichment focus on biological interpretation downstream of genetic analysis. We focus on the genetic analysis itself.
- *Spatial transcriptomics*: Unless directly connected to QTL analysis, spatial methods are beyond our scope.
- *Rare variant analysis*: Burden tests, SKAT, and related methods are well-motivated genetics problems but typically either lack rigorous statistical properties or are methodologically similar to common variant methods without new generative modeling insights. For theoretical discussion of rare vs. common variant disease etiology, readers should refer to population genetics work from Jonathan Pritchard and Guy Sella.
