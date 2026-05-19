# Datasets

This document describes all datasets used to evaluate the Empirical Bayes Calibration (EBC) check. The calibration procedure generates pseudo-data from the fitted empirical Bayes model and compares it to the observed data via a goodness-of-fit test using the energy statistic, yielding a Monte Carlo p-value. A small p-value (< 0.05) indicates that the fitted model is miscalibrated.

A summary of all results is provided in the table below, followed by a description of each dataset.



## Dataset Descriptions

### 2.2 Swirl

Source: Smyth et al. (2004), Example 9.1. Distributed as part of the `marray` R package.

Description: The goal of the Swirl experiment is to identify genes with altered expressions in the Swirl mutant compared to wild-type zebrafish. The observed data are pairs $(X_i, S_i^2)$ for $n = 8{,}448$ genes, where $X_i$ is the contrast for gene $i$ and $S_i^2$ is its standard error. Hyperparameters were estimated using the LIMMA package following Smyth's empirical Bayes model.

---

### 2.3 ApoAI

Source: Smyth et al. (2004), Example 9.2.

Description: This experiement compared 8 ApoAI knockout mice with 8 control mice by obtaining the target mRNA from the liver tissue from each of these 16 mice. The goal of this experiment is to determine how ApoAI defiiency affects the action of other genes in the liver. The gene expression microarray dataset has $n = 6{,}384$ genes. The same pre-processing pipeline and empirical Bayes procedure as the Swirl dataset were applied.

---

### 2.4 Chemical Abundance Ratios

Source: Soloff et al. (2024), Example 5.2.

Description: Chemical abundance ratios for the red clump (RC) stars given in the DR14 APOGEE red clump catalogue were examined. Following the same preprocessing steps as Ratcliffe et al. (2020), we have bivariate observations $X_i = (\text{Mg/Fe}_i, \text{Si/Fe}_i) \in \mathbb{R}^2$ and known heteroskedastic measurement error covariances $\Sigma_i$. The dataset contains n=29500 observations.

---

### 2.5 Math Scores in US Public Schools

Source: Soloff et al. (2024), Example 5.3.

Description: The Education Longitudinal Study (ELS) of 2002 contains math test scores and normalized socio-economic status (SES) of 10th grade students. The survey has n=100 different large public high schools on a total of $\sum_{i=1}^{n}N_i = 1993$ children. A d=2 hierarchical regression setting where $y_{ij}$ represents the math score of student j in school i, and $X_{ij}=[1,SES_{ij}]$ contains the corresponding SES score and an intercept term. The number of stduents $N_i$ surveyed in each school greatly varies, ranging from 4 to 32 with a median of 20 students.

---

### 2.5 Bball (Baseball)

Source: Gu et al 2017.

Description: Collected from ESPN, the data contains average monthly number of at bats and hits for all U.S Major League Baseball players from the regular seasons of 2002-2011, as well as an indicator of whether the player is a pitcher. We focus only on non-pitchers in the data analysis, leaving us with $n = 898$ MLB players and 9,199 observations. The batting averages are transformed as: $y_{it} = \arcsin\sqrt{(H_{it} + 0.25)/(N_{it} + 0.5)}$ where $H_{it}$ denotes the number of hits player $i$ achieves in period $t$, and $N_{it}$ denotes the number of at bats of the player in the same period. The $y_{it}$'s are assumed to be Gaussian with means $\mu_i = arcsin(\sqrt{p_i})$, where $p_i$ is the individual specific batting success probability and variances $\theta_i v_{it}^2 = \theta_i / (4N_{it})$.

---

### 2.6 Norberg

Source: REBayes R package.

Description: This dataset is available from the REBayes package as data("Norberg"). It consists of a portfolio of Norwegian workmen's gorup life insurance policies with $n = 72$ occupational groups. For each group $i$, $X_i$ denotes the number of deaths and $E_i$ is the corresponding exposure (number of years exposed to risk). The model assumes $X_i \mid \theta_i \sim \text{Poisson}(\theta_i E_i)$ with an unknown mixing distribution $F$ estimated via NPMLE. Observations are embedded as $Z_i = (\log E_i, \log(1 + X_i))$ for the energy distance computation.

---

### 2.7 Flies

Source: REBayes R package.

Description: This dataset is available from the REBayes package as data("flies"). It is a survival dataset with $n = 19,072$ observed lifetimes, with individual lifetimes $T_i > 0$ being modeled with a Weibull frailty mixture: $T_i \mid \theta_i \sim \text{Weibull}(\alpha, \theta_i)$. The frailty parameters $\theta_i \sim G$ drawn from an unknown mixing distribution estimated via the Kiefer-Wolfowitz NPMLE.

---

### 3. Chen 2025 — Income Data (Black and White)

Source: Chen (2025).

Description: We have Census tract-level estimates of mean income rank for Black children. Each observation $(Y_i, \sigma_i)$ consists of an estimated income rank and its standard error for census tract $i$. The empirical Bayes model assumes $Y_i \mid \theta_i, \sigma_i \sim N(\theta_i, \sigma_i^2)$, where $\theta_i = m_0(\sigma_i) + s_0(\sigma_i)\tau_i$ depends on precision $\sigma_i$ through unknown functions $m_0(\cdot)$ and $s_0(\cdot)$. The sample sizes is $n = 13{,}630$. 

---

### 4. Stars — Denoising / Jaffe

Source: Jaffe et al. (astronomy dataset).

Description: $n=2000$ stars were sampled from an astronomical dataset with $2.7 x 10^4$ stars. Latent signals $\Theta_i \in \mathbb{R}^d$ are observed with additive Gaussian noise $X_i \mid \Theta_i \sim N(\Theta_i, \Sigma_i)$, where $\Sigma_i$ are known diagonal measurement error covariance matrices. The unknown prior $G$ is assumed to be a Gaussian location mixture model and is estimated via heteroskedastic multivariate NPMLE. Additionally, a variance-constrained empirical Bayes denoiser is applied to correct for variance shrinkage.

---

### 6.1 Shakespeare (Deconvolver)

Source: Shakespeare's complete canon.

Description: Word frequency counts from Shakespeare's collected works are taken. $X_i$ denotes the total number of times distinct word $i$ appears in the canon, modeled as $X_i \mid \Theta_i \sim \text{Poisson}(\Theta_i)$. The observed data used for fitting are frequency counts $y_k = $ # $\{i : X_i = k\}$ for $k = 1, \ldots, 100$, giving $n = 30{,}688$ distinct words. The prior is estimated via g-modeling on a discrete log-grid using a spline-based exponential family.

---

### 6.2 Surgery (Deconvolver)

Source: Intestinal cancer surgery records.

Description: The data has $n = 844$ cancer patients. For each patient $i$, the observed pair $(n_i, X_i)$ records the number of satellite nodes removed ($n_i$) and the number found to be malignant ($X_i$). The model assumes $X_i \mid \theta_i \sim \text{Binomial}(n_i, \theta_i)$, where $\theta_i$ is the patient-specific probability of malignancy, and $\theta_i \sim G$ is estimated empirically.

---

### 7. Statin46 (EBGG / General-Gamma)

Source: Pharmacovigilance adverse event reporting system.

Description: The data is an adverse event–drug contingency table for six statin drugs, with $n = 276$ adverse event–drug pairs. For each pair $(i, j)$, $N_{ij} \sim \text{Poisson}(E_{ij} \lambda_{ij})$, where $E_{ij}$ is a null expected count and $\lambda_{ij}$ is a latent signal parameter estimated via a flexible gamma-mixture empirical Bayes prior. Four different EB methods are compared: GPS, KM, Efron, and General-Gamma.

---

### 9. Hockey (Mindist)

Source: 2017–18 NHL season goal records.

Description: Goal totals $Y_i$ for $n = 818$ NHL players in the 2017–18 season, modeled as $Y_i \mid \theta_i \sim \text{Poisson}(\theta_i)$ with an unknown prior $G$ on player scoring rates. Three nonparametric minimum-distance estimators of $G$ are compared: NPMLE (KL divergence), Hellinger distance, and $\chi^2$ divergence. The fitted EB rules are used to predict 2018–19 goal totals.

---

### 11. Species Richness Datasets

Source: Various ecological surveys.

Description: Five ecological count datasets where $Y_i$ denotes the number of times species $i$ is detected across $s$ sampling efforts. Each is modeled as $Y_i \mid \lambda_i \sim \text{Poisson}(\lambda_i)$ with the prior $G$ estimated via a nonparametric g-modeling approach. The calibration test compares observed positive-count samples to pseudo-samples from the fitted model.

| Dataset | n | P-value |
|---------|---|---------|
| Traffic | 1,621 | 0.768 |
| Microbial | 514 | 0.668 |
| Brain Vessel | 365 | 0.712 |
| Hard Candy | 354 | 0.730 |
| Death Notices | 934 | 0.728 |

---

### 12. High-Dimensional Microarray (NPMLE)

Source: Golub et al. and Lung cancer gene expression studies.

Description: Two high-dimensional gene expression datasets with binary class labels. For each gene $j$ and class $k \in \{0, 1\}$, the class-specific sample mean $\bar{X}_j^k \sim N(\mu_j^k, 1/n_k)$ is modeled with an unknown prior $G_k$ estimated separately per class via the Kiefer–Wolfowitz NPMLE.

| Dataset | n (genes) | P-value |
|---------|-----------|---------|
| Golub | 7,129 | 0.029, 0.057 |
| Lung | 12,533 | 0.079, 0.086 |

---

### 13. Binary Car (NPMLE)

Source: Modal choice transportation survey.

Description: A binary discrete-choice dataset with $n = 842$ observations. The binary outcome $y_i \in \{0, 1\}$ encodes a modal choice, modeled via a random coefficients probit: $y_i = \mathbf{1}\{x_i^\top \beta_i > 0\}$ with latent random coefficients $\eta_i \sim F$. The unknown mixing distribution $F$ is estimated by NPMLE.

---

### 14. NPMLE Baseball

Source: MLB first-half batting records.

Description: First-half at-bats $A_j$ and hits $H_j$ for $n = 565$ MLB players. The hierarchical model assumes $A_j \mid (\lambda_j, \pi_j) \sim \text{Poisson}(\lambda_j)$ and $H_j \mid (A_j, \lambda_j, \pi_j) \sim \text{Binomial}(A_j, \pi_j)$, with the joint latent parameter $(\lambda_j, \pi_j) \sim G_0$ estimated via an approximate bivariate NPMLE on a fixed grid.

---

### 15. Breast Microarray (NPMLE Mixture Joint)

Source: Breast cancer gene expression study.

Description: Gene expression data with a binary outcome (response/no response) for $p = 22{,}283$ genes across $n$ subjects. For each gene $j$, the pair of class-specific sample means $Z_j = (\bar{X}_{j0}, \bar{X}_{j1})$ is modeled jointly as $Z_j \mid \Theta_j \sim N_2(\Theta_j, \text{diag}(1/n_0, 1/n_1))$, with a bivariate prior $G$ estimated via NPMLE on a finite grid.

---

### 16. STAR (Gauss and NPMLE)

Source: Tennessee Student/Teacher Achievement Ratio (STAR) class-size experiment.

Description: School-level estimated treatment effects $\hat{\theta}_{ij}$ across $J = 4$ strata (class-size conditions), with estimated variances $S_{ij}^2$, for $n = 5{,}771$ schools. The model assumes $\hat{\theta}_i \mid \theta_i \sim N(\theta_i, \Sigma_i)$ with the prior $G$ estimated via either a parametric Gaussian model or NPMLE.

---

### 17. Motorcycle (SMASH)

Source: Motorcycle crash simulation experiment.

Description: A simulated motorcycle crash acceleration signal with $T = 128$ time points after preprocessing (median aggregation, ordering, and padding to a power of 2). The heteroskedastic Gaussian model $y_t \mid \mu_t, \sigma_t^2 \sim N(\mu_t, \sigma_t^2)$ is fitted using SMASH, an empirical Bayes wavelet shrinkage method that estimates both the mean signal and the variance function iteratively.

---

### 18. wOBA Baseball (ebnm)

Source: MLB player performance records.

Description: Weighted on-base average (wOBA) observations $(x_i, s_i)$ for $n = 688$ MLB players, where $x_i$ is the observed wOBA and $s_i$ is its estimated standard error (treated as known). The normal means model $x_i \mid \theta_i \sim N(\theta_i, s_i^2)$ is fitted with three different empirical Bayes prior families (e.g., normal, unimodal, nonparametric mixture) via the `ebnm` package.

---

### 19. neuralG — Arizona and FHS

Source: Arizona population study (Poisson); Framingham Heart Study (FHS, measurement error model).

Description: Two datasets used to evaluate the `neuralG` g-modeling empirical Bayes approach, each with $n = 1{,}798$ observations.

- Arizona: Count data modeled as $Y_i \mid \lambda_i \sim \text{Poisson}(\lambda_i)$ with an unknown prior $G$ estimated via neural-g on a discrete grid.
- FHS: Continuous measurements with replication noise, modeled as $Y_i \mid \theta_i \sim N(\theta_i, \sigma^2_{\text{obs}})$ (two replicates averaged), with $G$ estimated similarly.
