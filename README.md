# MixIt
The R source code for [Jang, Jaehyuk, Suehyun Kim, and Kwonsang Lee. "Improving Causal Estimation by Mixing Samples to Address Weak Overlap in Observational Studies." *arXiv preprint arXiv:2411.10801 (2024)*.](https://arxiv.org/abs/2411.10801)

### Abstract
In observational studies, the assumption of sufficient overlap (positivity) is fundamental for the identification and estimation of causal effects. Failing to account for this assumption yields inaccurate and potentially infeasible estimators. To address this issue, we introduce a simple yet novel approach, \textit{mixing}, which mitigates overlap violations by constructing a synthetic treated group that combines treated and control units. Our strategy contributes to weighting literature with several advantages. First, it improves the accuracy of the estimator by preserving unbiasedness while reducing variance. This phenomenon results from the shrinkage of propensity scores in the mixed sample, which enhances robustness to poor overlap, but remains effective regardless of the overlap level. Second, it enables direct estimation of the target estimand without discarding extreme observations or modifying the target population. Third, the mixing approach is highly adaptable to various weighting schemes, including contemporary methods such as entropy balancing. Fourth, it serves as an indirect method to diagnose model misspecification. We illustrate its empirical performance through extensive simulation studies. We also introduce a practical guidance for diagnosing and tackling limited overlap with an analysis of chronic obstructive pulmonary disease (COPD) dataset provided by Samsung Medical Center. Through our strategy, we investigate whether small airway disease (SAD) precedes emphysema progression among early-stage COPD patients.

### Files Description
| File Name | Description |
| --------- | ----------- |
| functions.R | The core functions to implement mixing strategy |
| utils.R | The utility functions to aid convenience in simulation studies |
| sim1.R | The simulation study of section 4.1 in the main paper |
| sim2.R | The simulation study of section 4.2 in the main paper |
| copd.R | The real data analysis of section 5 in the main paper |
| adv_mix.R | The simulation study of section 6.1 in the main paper |
| data.R | The COPD dataset from Samsung Medical Center |
| smc.boot.mat.RData | The bootstrap estimates of the COPD study |

### Simple Mixing Strategy
The *simple mixing strategy* estimates the treatment effect of the treated (ATT) through  the mixed samples drawn from $h^*$, not from the observed samples drawn from the density $h$. This strategy is developed in order to mitigate the weak overlap issue in causal inference, particularly in weighting literature. For $\delta \in (0,1)$, we define the simple mixed distribution as

$$
\begin{aligned}
    h^*_1 &= (1-\delta) h_1 + \delta h_0 \\
    h^*_0 &= h_0
\end{aligned}
$$

This implies that the synthetic treated group is a mixture of the original treated and control group with ratio $1-\delta: \delta$. On the other hand, the synthetic control group is identical to the original. Leveraging the mixed distribution, we can estimate the ATT with two different algorithms: M-estimation-based and mixing-algorithm-based, each having its own advantages and disadvantages. For further details, see the main paper.

The following `mixit()` function takes argument of either approaches to estimate the ATT.

```{R}
mixit(
  formula,            # treatment ~ covariates
  outcome,            # character: name of the observed-outcome column
  data,               # data.frame with outcome, treatment, covariates
  delta,              # mixing proportion (0 < delta < 1)
  by,                 # "M-estimation"  or  "Algorithm"
  weight  = NULL,     # if by = "Algorithm": "glm", "ebal", or "CBPS"
  M       = 200,      # Monte-Carlo replicates for the mixing algorithm
  ...                 # passed to  weightit()  or  CBPS()
)
```

| Argument  | Purpose                                                                              | Details / default                                                                                                               |
| --------- | ------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------- |
| `formula` | Defines the treatment indicator and baseline covariates.                             | Standard R formula, e.g. `treated ~ age + bmi + smk + Emph + interval`.                                                                           |
| `outcome` | Column name of the observed outcome `Y`.                                             | Character scalar.                                                                                                               |
| `data`    | Original sample. Must contain `outcome`, treatment, and all covariates.              | `data.frame` or tibble.                                                                                                         |
| `delta`   | Mixing proportion $\delta\in(0,1)$.                                                  | $\delta=0$ ⇒ no mixing; larger $\delta$ ⇒ more controls mixed into the treated group.                                           |
| `by`      | Implementation route.                                                                | `"M-estimation"` or `"Algorithm"`                                |
| `weight`  | Weighting method for `"Algorithm"`.                                                  | `"glm"` (propensity GLM), `"ebal"` (entropy balancing), or `"CBPS"` (covariate-balancing PS). Ignored if `by = "M-estimation"`. |
| `M`       | Number of resampled mixed datasets when `by = "Algorithm"`.                          | Integer ≥ 1, default `200`.                                                                                                     |
| `...`     | Extra arguments forwarded to `weightit()` (for `"ebal"`) or `CBPS()` (for `"CBPS"`). |                                                                                                                                 |


#### M-estimation-based approach (`by = "M-estimation"`)
This approach is currently only available for designing the original propensity scores with logistic model.

```{R}
source('functions.R')      # load functions
load(file = 'data.RData')  # load some data

# --- specify inputs ----
form   <- treated ~ age + bmi + smk + Emph + interval
d  <- 0.10

# M-estimation version --------------------------------------
mipw_M <- mixit(
  formula = form,
  outcome = "Emph_outcome",
  data    = data,
  delta   = d,
  by      = "M-estimation"
)

mipw_M$est      # point estimate
sqrt(mipw_M$V)  # sandwich SE
mipw_M$weights  # weights
```



#### Mixing-algorithm-based approach (`by = "Algorithm"`)
This approach currently extends the mixing to be implemented onto `glm`, `ebal` and `CBPS` weighting schemes.

```{R}
source('functions.R')      # load functions
load(file = 'data.RData')  # load some data

# --- specify inputs ----
form   <- treated ~ age + bmi + smk + Emph + interval
d  <- 0.10

# Mixing algorithm + entropy balancing ----------------------
meb_alg <- mixit(
  formula = form,
  outcome = "Emph_outcome",
  data    = data,
  delta   = d,
  by      = "Algorithm",
  weight  = "ebal",
  M       = 200              # more replicates for stability
)

meb_alg$est             # point estimate
head(meb_alg$weights)   # inspection of averaged weights
```
