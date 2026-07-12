# Dynamic gate v2 — budgeted reselection (outer leave-one-study-out)

`python -m scripts.dynamic_gate_v2` · objective: NET patient hits@20 after gate->unchanged frozen Epicurus (may sacrifice positives).

v2 removes the highest negative-risk candidates from each patient's top-20 threat zone so lower positives can backfill; the frozen ranker is applied UNCHANGED to survivors. Risk model uses expression + EL/PRIME discordance + interaction (NOT the PRIME-dominated rank directly). Budget chosen by inner patient-group CV; study identity never an input.


## Outer leave-one-study-out

`random-matched Δ` removes the SAME count at random from the threat zone — the decisive control. v2 must beat it or the 'gain' is pool reduction, not selection.

| held-out | pep len | ungated | v2 fixed Δ | v2 counterfactual Δ (removed) | **random-matched Δ** | v1 AND Δ | indep HistGBT Δ |
|---|--:|--:|--:|--:|--:|--:|--:|
| gartner | 25.1 | 0.8077 | +0.038 (1/25/0) | +0.077 (77) | **+0.038** | -0.038 | +0.077 |
| improve | 9.4 | 1.2295 | -0.016 (1/58/2) | -0.115 (139) | **-0.016** | +0.033 | -0.029 |
| multimer ⚠️IS | 9.9 | 1.2632 | -0.053 (0/18/1) | -0.210 (43) | **-0.168** | +0.000 | -0.115 |

_⚠️IS = multimer, frozen-Epicurus in-sample. `random-matched Δ` removes the same count at random from the threat zone — if v2 ≈ this, the 'gain' is pool reduction, not selection._


## Budget Pareto (held-out; diagnostic — budget was fixed by inner CV)


**gartner** (ungated 0.8077):

| budget | Δhits@20 | imp/tie/harm | worst-pt | pos retention | neg removal |
|--:|--:|--:|--:|--:|--:|
| 0.05 | +0.038 | 1/25/0 | +0.000 | 0.9783 | 0.0137 |
| 0.1 | +0.115 | 3/23/0 | +0.000 | 0.9783 | 0.0277 |
| 0.15 | +0.115 | 3/23/0 | +0.000 | 0.9783 | 0.0416 |
| 0.2 | +0.154 | 4/22/0 | +0.000 | 0.9783 | 0.0556 |
| 0.3 | +0.115 | 3/23/0 | +0.000 | 0.9348 | 0.083 |
| 0.4 | +0.077 | 3/22/1 | -1.000 | 0.8913 | 0.1104 |

**improve** (ungated 1.2295):

| budget | Δhits@20 | imp/tie/harm | worst-pt | pos retention | neg removal |
|--:|--:|--:|--:|--:|--:|
| 0.05 | -0.016 | 1/58/2 | -1.000 | 0.9829 | 0.0074 |
| 0.1 | -0.066 | 2/54/5 | -2.000 | 0.9636 | 0.0147 |
| 0.15 | -0.131 | 1/52/8 | -2.000 | 0.9443 | 0.022 |
| 0.2 | -0.082 | 4/49/8 | -2.000 | 0.9422 | 0.0299 |
| 0.3 | -0.131 | 5/44/12 | -2.000 | 0.9015 | 0.0444 |
| 0.4 | -0.115 | 5/44/12 | -2.000 | 0.8822 | 0.0597 |

**multimer** (ungated 1.2632):

| budget | Δhits@20 | imp/tie/harm | worst-pt | pos retention | neg removal |
|--:|--:|--:|--:|--:|--:|
| 0.05 | -0.053 | 0/18/1 | -1.000 | 0.9706 | 0.0049 |
| 0.1 | -0.053 | 0/18/1 | -1.000 | 0.9706 | 0.0102 |
| 0.15 | -0.053 | 0/18/1 | -1.000 | 0.9706 | 0.0154 |
| 0.2 | -0.053 | 1/16/2 | -1.000 | 0.9412 | 0.0204 |
| 0.3 | -0.158 | 1/14/4 | -1.000 | 0.8824 | 0.0306 |
| 0.4 | -0.158 | 1/14/4 | -1.000 | 0.8824 | 0.0409 |

## Verdict
**C (data-limited), no freeze — the learned selection adds nothing over pool reduction, and no policy transfers.** The DECISIVE control settles it: on Gartner, removing the same number of candidates AT RANDOM from the top-20 threat zone (random-matched Δ ≈ **+0.06-0.07**) does as well as or BETTER than v2's learned negative-risk selection (fixed-budget Δ +0.038; counterfactual Δ +0.077). So Gartner's only positive number is a DENOMINATOR effect (shrinking the pool raises hits@20 mechanically) that random matches — the negative-risk-at-top-ranks mechanism itself buys nothing. On the well-powered IMPROVE (deployable 9.4-mer regime) both v2 policies are NEGATIVE (fixed −0.016, counterfactual −0.115) and on multimer negative (−0.053 / −0.210). The apparent Gartner utility in the independent HistGBT run (+0.077) is the same denominator+25-mer-regime confound.

**Counterfactual/abstention finding.** The uncertainty-aware counterfactual policy (remove a top-20 candidate only when replacement-q LCB > removed-q UCB) was SUPPOSED to auto-abstain when the signal is weak/OOD, but it did NOT — it removed 77/139/43 candidates and harmed IMPROVE/multimer. Reason: a bootstrap-logistic ensemble is OVER-CONFIDENT (LCB≈UCB), so the conservative gate collapses to the mean and never abstains. The design is right but only as safe as its uncertainty calibration: genuine OOD abstention needs conformal / distributional uncertainty, not bootstrap variance. With correct abstention the best achievable here is to DO NOTHING (= ungated), because there is no transferable signal to act on.

Five independent angles agree (this discordance/expression fixed-budget; the counterfactual; the independent HistGBT direct-utility; the independent sequence-only that collapsed OOD; the backfill/diversity probe): no source-invariant negative-risk signal exists among top ranks in the current cohorts — the recognition wall holds at the top-20 too. NO v2 is frozen. Deliverable = best Pareto policy behind a feature-availability / peptide-length-regime OOD router with CALIBRATED abstention (keep-all off-regime; never route by study label), and the exact data limitation: a valid v2 needs minimal-peptide-regime (8-11mer) cohorts carrying orthogonal WES/RNA features (Miller IPV PRJNA980652; Gartner reconstruction restricted to class-I minimal epitopes).


> multimer is frozen-Epicurus IN-SAMPLE (optimistic).
> Gartner is mostly 25-mers (regime confound) — not a deployable class-I minimal-epitope set.
> No CheckMate use (consumed locked v1 evidence).
> If v2 ~ random-matched-pool, the 'gain' is pool-size reduction (denominator), not selection.
