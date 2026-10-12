# AI/deep-learning benchmark rationale

A focused literature check identified Cox-compatible deep survival models, interpretable pathway-based networks, two-stage Cox-nnet designs and newer patient-similarity graph neural networks as the most relevant AI directions for transcriptomic prognosis work.

The implemented benchmark therefore used a small DeepSurv/Cox-nnet-style multilayer perceptron trained with the negative Cox partial log-likelihood. To preserve the deployability-audit framing, the feature space was restricted to the published-signature gene universe rather than unrestricted whole-transcriptome discovery.

Key sources used for design context:
- TRIPOD+AI statement, BMJ 2024, doi:10.1136/bmj-2023-078378.
- PROBAST+AI, BMJ 2025, doi:10.1136/bmj-2024-082505.
- Huang et al., BMC Medical Genomics 2020, doi:10.1186/s12920-020-0686-1.
- Hao et al., BMC Medical Genomics 2019, doi:10.1186/s12920-019-0624-2.
- Zhan et al., NAR Genomics and Bioinformatics 2021, doi:10.1093/nargab/lqab015.
- Mao et al., Briefings in Bioinformatics 2025, doi:10.1093/bib/bbaf108.

Benchmark result: DeepSurv external median C-index 0.476; linear Cox external median C-index 0.544.

Interpretation: This is an exploratory AI augmentation benchmark, not a newly deployable clinical prediction model.