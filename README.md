# SEAGLE

The explosion of biobank data offers immediate opportunities for gene-environment (GxE) interaction studies of complex diseases because of the large sample sizes and rich collection in genetic and non-genetic information. 
However, the extremely large sample size also introduces new computational challenges in GxE assessment, especially for set-based GxE variance component (VC) tests, a widely used strategy to boost overall GxE signals and to evaluate the joint GxE effect of multiple variants from a biologically meaningful unit (e.g., gene). 

We present SEAGLE, a Scalable Exact AlGorithm for Large-scale Set-based GxE tests, to permit GxE VC test scalable to biobank data. SEAGLE employs modern matrix computations to achieve the same “exact” results as the original GxE VC tests, and does not impose additional assumptions nor relies on approximations. SEAGLE can easily accommodate sample sizes in the order of $10^5$, is implementable on standard laptops, and does not require specialized equipment. 

## Installation

To install the latest stable version from CRAN:

  ```{r}
install.packages('SEAGLE')
```

To install the latest development version from GitHub:

  ```{r}
# install.packages("devtools")
devtools::install_github('jocelynchi/SEAGLE')
```

## Getting Started

We've included four examples on how to use the `SEAGLE` software.  

1. The [first example](https://jocelynchi.github.io/SEAGLE/articles/example1.html) shows how to use `SEAGLE` when the user inputs phenotype, covariate, and genotype data from .txt files.  
2. The [second example](https://jocelynchi.github.io/SEAGLE/articles/example2.html) shows how to use `SEAGLE` when the user has a genetic marker matrix **G** that is already in matrix form.  
3. The [third example](https://jocelynchi.github.io/SEAGLE/articles/example3.html) shows how to use `SEAGLE` when the user has genotype data from GWAS or next generation sequencing studies.
4. The [fourth example](https://jocelynchi.github.io/SEAGLE/articles/example4.html) shows how to use `SEAGLE` for chromosome-wide analysis when the user has genotype data from GWAS studies.

## Citing SEAGLE

The accompanying journal manuscript for `SEAGLE` can be found at [arXiv:2105.03228](https://arxiv.org/abs/2105.03228).  To cite the `SEAGLE` software, please use the following BibTeX entry.

```
@misc{seagle,
  author = {Jocelyn T. Chi and Ilse C. F. Ipsen and Tzu-Hung Hsiao and Ching-Heng Lin and Li-San Wang and Wan-Ping Lee and Tzu-Pin Lu and Jung-Ying Tzeng},
  title = {SEAGLE: A Scalable Exact Algorithm for Large-Scale Set-based GxE Tests in Biobank Data},
  year = {2021+},
  arxiv = {https://arxiv.org/abs/2105.03228},
  howpublished = {arXiv:2105.03228 [stat.CO]}
}
```

## Acknowledgments

Many thanks to Yueyang Huang for his help with generating the example data and PLINK1.9 code for the tutorials.

