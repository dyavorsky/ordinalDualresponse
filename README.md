# ordinalDR

Ordinal dual-response choice models for conjoint analysis.

Choice-based conjoint studies increasingly pair each forced choice with an
**ordinal** purchase-likelihood follow-up ("how likely are you to buy the
option you selected?" on a 5- or 7-point scale) instead of the binary
buy/no-buy question of the dual-response literature. `ordinalDR` provides the
statistical model for that design.

## The model in five equations

With inside-good utilities $u_{j} = x_j'\beta + \eta_j$, $\eta_j \sim$ i.i.d.
Gumbel(0,1), and an outside-good shock $\eta_0$ the consumer has not resolved:

1. **Forced choice** $\;j^* = \arg\max_j u_j$: standard MNL over inside goods,
   $\Pr(j) = e^{V_j} / S$, $S = \sum_k e^{V_k}$.
2. **Max-stability**: $u^* = \max_j u_j \sim$ Gumbel$(\bar\mu, 1)$ with
   $\bar\mu = \ln S$, independent of *which* good attained it.
3. **Ordinal report**: cumulative-link model
   $\Pr(y \le w) = G(c_w - \bar\mu)$ — Gumbel link $G$ if the consumer
   reports her purchase *probability* in bins (Model A), logistic link if she
   grades the resolved comparison $u^* - \eta_0$ (Model B).
4. **Joint likelihood** = MNL $\times$ interval probability (a consequence of
   max-stability, not an independence assumption).
5. **Nesting**: at $W = 2$, Model B collapses *exactly* to the binary
   dual-response "Unified Model" of Hung et al. (2025) / DR-AnyMax of Diener
   et al. (2006) — equivalently, standard MNL with a no-choice constant.

## What's in the package

- `dr_design()`, `dr_panel_design()`, `dr_simulate()` — designs and
  behavioral-process simulation (heterogeneous preferences and cut points,
  DR-2Max and slope alternatives, non-Gumbel shocks, task-order drift)
- `dr_mle()` — aggregate maximum likelihood with tail-stable link math
- `dr_hb()` — hierarchical Bayes by RW-Metropolis-within-Gibbs, with common
  or respondent-specific cut points; `dr_hb_betai()`, `dr_holdout_ll()`
- `dr_fisher()`, `dr_omega()`, `dr_optimal_cuts()`, `dr_equal_share_cuts()`,
  `dr_power()` — information and scale-design tools
- `dr_spec_tests()` — the unit-slope ($\lambda = 1$) and conditioning
  ($\delta = 0$ vs DR-2Max's $\delta = 1$) specification tests
- `dr_lmd_nr()`, `dr_bridge_lmd()`, `dr_laplace_lmd()`, `dr_waic()` — model
  comparison
- `dr_demand_dual()`, `dr_demand_constant()`, `dr_nash_prices()` —
  counterfactual pricing

## Installation

```r
# development version (private during paper review)
remotes::install_github("dyavorsky/ordinalDR")
```

## Quick start

```r
library(ordinalDR)
beta <- c(a2 = 0.8, a3 = -0.5, b2 = 0.4, b3 = 1.0, price = -0.9)
cut  <- c(-1.0, 0.2, 1.2, 2.2)                  # 5-point scale
des  <- dr_design(3600, J = 4, seed = 1)
dat  <- dr_simulate(des, matrix(beta, 1, 5), cut, model = "B", seed = 2)
dr_mle(dat)
dr_spec_tests(dat)
```

See the vignettes: `minimal-example`, `aggregate-study`, `hierarchical`.

## Reference

Bhalerao, Yavorsky & Zheng, "Outside Good Uncertainty: Ordinal Dual Response
in Choice-Based Conjoint Analysis" (working paper).
