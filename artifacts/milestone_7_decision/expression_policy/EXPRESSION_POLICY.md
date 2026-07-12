# RNA-expression ranking-policy analysis (label-blind, development cohorts)

> Policy `expression-policy-analysis-1.0.0`, k=20. Frozen decision: **confidence_only** — prime_only_incumbent (genuine PRIME; expression does NOT move rank).

> Development cohorts are interpreted WITHIN their own denominators and NEVER pooled. The protected incumbent is lossless generation + genuine PRIME.

## Structural: recognition vs within-patient expression quartile (per cohort)

| cohort | Q1lo | Q2 | Q3 | Q4hi | positives in low-expr half |
|---|---|---|---|---|---|
| cd8_multimer | 0.0029 | 0.003 | 0.0045 | 0.0063 | 12/34 |
| gartner_nci | 0.0 | 0.0006 | 0.0031 | 0.0175 | 1/46 |
| improve_srhgroup | 0.024 | 0.0254 | 0.0318 | 0.0252 | 213/467 |

## Per-cohort no-regression vs PRIME incumbent (recall@20; ✅ no-regression / ❌ regresses)

| cohort | prime_only | expr_rank_penalty | soft_saturating | portfolio_reserve |
|---|---|---|---|---|
| cd8_multimer | 0.6176 ✅ | 0.4706 ❌ | 0.6176 ✅ | 0.6176 ✅ |
| gartner_nci | 0.3696 ✅ | 0.5652 ✅ | 0.3696 ✅ | 0.3043 ❌ |
| improve_srhgroup | 0.1777 ✅ | 0.1606 ❌ | 0.1777 ✅ | 0.1777 ✅ |

## Decision

- **Chosen: confidence_only** — prime_only_incumbent (genuine PRIME; expression does NOT move rank)
- Expression role: confidence annotation (stratum label) + optional off-by-default route-dependent portfolio reserve for reachability
- No-regression everywhere: prime_only_incumbent=YES, expr_rank_penalty=no, soft_saturating=YES, portfolio_reserve=no

Per development cohort (never pooled): expr_rank_penalty regresses [cd8_multimer, improve_srhgroup] while only helping [gartner_nci]; portfolio_reserve regresses [gartner_nci] within a fixed top-20 budget. The two no-regression-everywhere forms are prime_only (confidence-only) and soft_saturating, and on these cohorts soft_saturating is IDENTICAL to prime_only (it only demotes candidates that are BOTH low-presentation and low-expression, which never occupy the top-20) — i.e. letting expression move rank buys NO measurable benefit and risks harm. Decision: expression is CONFIDENCE-ONLY in the score; lossless+PRIME stays the protected incumbent ranking. The soft-saturating guard is frozen as an equivalent, route-dependent no-op-on-top-20 safeguard, and the portfolio reserve is retained as an OPTIONAL off-by-default reachability tool (it trades displacement for stratum coverage). Because PRIME already keeps high-presentation candidates irrespective of expression, this preserves reachability of low-expression recognized candidates.

## Sid descriptive (post-freeze, n=3 — NOT a gate)

- frozen confidence-only (= lossless+PRIME): recall@20 3/3 (ASPM, DYNC1H1, MAP2)
- counterfactual expr rank penalty: recall@20 2/3 (ASPM, DYNC1H1)
- counterfactual soft-saturating: recall@20 3/3 (ASPM, DYNC1H1, MAP2)

Confidence-only (= lossless+PRIME) keeps all 3 recognized mutations in the top-20; the expression rank penalty demotes low-expression MAP2. Consistent with the development no-regression finding. Descriptive on n=3; nothing tuned to these labels.

