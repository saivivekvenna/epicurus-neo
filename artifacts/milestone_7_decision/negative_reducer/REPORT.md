# Risk-controlled negative reducer — nested LOSO (non-Sid)

_DEV 122 patients / 29391 rows. NO Sid/Miller. HGB available True (diagnostic-only). git 445208880f; sklearn 1.9.0; scipy 1.18.0; reproducible=True._


## Outer leave-one-study-out

| held-out | chosen (model,C,m) | inner Δ | test retention | CP-LB (pow) | neg removed | Δhits@20 | matched-rand | beats | pos-rank↑ |
|---|---|--:|--:|--:|--:|--:|--:|:--:|--:|
| improve | NULL,None,m=0 | +0.000 | 1.0 | 0.9936 (y) | 0.0 | +0.000 | +0.000 | n | +0.00 |
| gartner | NULL,None,m=0 | +0.000 | 1.0 | 0.937 (n) | 0.0 | +0.000 | +0.000 | n | +0.00 |
| multimer | NULL,None,m=0 | +0.000 | 1.0 | 0.9157 (n) | 0.0 | +0.000 | +0.000 | n | +0.00 |

## Aggregate
Δhits@20 = **+0.0000** (matched-random +0.0000); bootstrap +0.0000 CI[0.0, 0.0] p>0=0.0; worst study +0.0000. Aggregate CP-LB 0.9945 (0/547 pos removed); IMPROVE CP-LB 0.9936 (0/467).


## §5 eligibility
- every_study_noncatastrophic: **True**
- aggregate_gain_beats_random: **False**
- any_negative_removal: **False**
- aggregate_cp_ge_0_95: **True**
- improve_cp_ge_0_95: **True**
- loso_eligible: **False**
- deploy_recipe_valid: **False**
- ELIGIBLE: **False**

## FROZEN: **NULL** — nested LOSO evidence did not pass CONTRACT §5 (see eligibility)
full-DEV selection {'model': 'NULL', 'C': None, 'm': 0, 'delta': 0.0}; SHA `6570612678b92e8a`. Honest negative.


**Verdict: NULL FROZEN (honest negative).** Sid/Miller locked; no application performed.
