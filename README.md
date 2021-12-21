# Boosted Dyadic Kernel Discriminant

## EECS6690-2021Fall Final Project

Reproduced and improved the Boosted Dyadic Kernel Discriminant. 

## Authors
* Yuqing Cao (yc3998)
* Fengyang Shang (fs2752)
* Chengbo Zang (cz2678)

## Description

Reproduced results of Boosted Dyadic Kernel Discriminant[1], introcuded bootstrap and cross validatioin, conducted experiments trying to improve the performance, impelmented and analyzed dimensionality reduction. 

## Getting Started

### Dependencies

A few of the key tools and environments are listed as follows. 
* RCurl
* pracma
* resample
* rdetools
* caTools

## Directory Orgnization

```
.
├── README.md     - Introduction
├── cv.R          - Crossvalidations & Parameter Selection
├── datasets.R    - Experiments on Defferent Datasets
├── dimred.R      - Dimesionality Reduction
├── hypercuts.R   - Model Definition
└── svm.R         - Comparisons With SVM

0 directories, 7 files
```

## References

[1] Moghaddam, Baback & Shakhnarovich, Gregory. (2003). Boosted Dyadic Kernel Discriminants.

[2] B. E. Boser, I. M. Guyon, and V. N. Vapnik. A training algorithm for optimal margin classifiers. In D. Haussler, editor, Proc. 5th Annual ACM Workshop on Computational Learning Theory, pages 144-152. ACM Press, 1992.

[3] C. L. Blake and C. J. Merz. UCI repository of machine learning databases. [https://www.ics.uci.edu/~mlearn/MLRepository.html], 1998.

[4] Gareth James, Daniela Witten, Trevor Hastie, Robert Tibshirani. An Introduction to Statistical Learning : with Applications in R. New York :Springer, 2013.

[5] Jolliffe Ian T. and Cadima Jorge 2016 Principal component analysis: a review and recent developments Phil. Trans. R. Soc. A.3742015020220150202.

[6] Course Material from CMU, CS Department, Computer Science Theory for the Information Age, Spring 2012.

