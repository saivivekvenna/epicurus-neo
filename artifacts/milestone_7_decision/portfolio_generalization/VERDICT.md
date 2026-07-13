# Portfolio generalization verdict

## Bottom line

The Hu_287 3/3 is real and reproducible under the frozen policy, but the current
two-patient evidence does **not** establish that the selector generally improves
recognition. It establishes a narrower, mechanistically coherent result:

> Epicurus v0.1 had useful **mutation-level ordering** on Hu_287, but ordinary
> peptide-level top-k selection hid it under many redundant peptide×HLA routes.
> The frozen mutation-diversity cap exposed that signal and produced 3/3.

On Sid, the selector diversifies the list but cannot repair weak mutation-level
ordering, so it does not change hits@20.

## Primary crossed result

| Patient | PRIME plain | PRIME + frozen selector | Epicurus plain | Epicurus + frozen selector |
|---|---:|---:|---:|---:|
| Hu_287 | 0/3 | 2/3 | 0/3 | **3/3** |
| Sid | **2/3** | **2/3** | 1/3 | 1/3 |

The route reserves contribute zero additional hits in both patients. The effect
comes entirely from the frozen per-mutation diversity cap.

## Why Hu_287 improves

Best-route ranks after collapsing each mutation to its best candidate:

| Recognized mutation | PRIME mutation rank | Epicurus mutation rank |
|---|---:|---:|
| CIITA `16:10907146:C:T` | 13 | **4** |
| TP53 `17:7673535:C:G` | 10 | **8** |
| SEC23B `20:18511052:C:T` | 11 | **5** |

Epicurus therefore put all three positive mutations inside its best eight
mutations. Its ordinary peptide top 20 nevertheless represented only **three
unique mutations** and contained **17 duplicate-route slots**. With at most two
routes per mutation, its selected set represented ten mutations and recovered
all three positives.

PRIME's mutation ordering was weaker (10th, 11th, and 13th), so the same cap-2
selector reached two positives. This is why the full 3/3 is an interaction:
Epicurus supplied better mutation-level ordering on this patient, and portfolio
selection stopped route redundancy from concealing it.

## Why Sid does not improve

| Recognized mutation | PRIME mutation rank | Epicurus mutation rank |
|---|---:|---:|
| DYNC1H1 | **3** | **1** |
| MAP2 | **9** | 71 |
| ASPM | 41 | 44 |

PRIME already selected DYNC1H1 and MAP2. Its ordinary top 20 represented 17
unique mutations, leaving little duplicate burden for the cap to correct.
Epicurus's selector increased unique-mutation coverage from 8 to 13, but MAP2
and ASPM were genuinely ranked too low at mutation level. Diversification cannot
turn the 44th- and 71st-ranked mutations into a defensible 20-slot choice.

## Sensitivity results

- At `k=10`, Hu_287 Epicurus improves 0/3 → 2/3; PRIME remains 0/3.
- At `k=30`, Hu_287 Epicurus improves 1/3 → 3/3 and PRIME 2/3 → 3/3.
- Sid remains unchanged by the selector at `k=10`, `20`, and `30`.
- A post-hoc cap of one route per mutation gives both scorers 3/3 on Hu_287 but
  still leaves Sid at PRIME 2/3 and Epicurus 1/3. It must not replace the frozen
  cap-2 headline.
- Caps 3, 5, or none weaken the Hu_287 effect, confirming that duplicate-route
  concentration—not route reserves—caused the gain.

## What is supported now

1. Treating a vaccine shortlist as a set-selection problem can materially change
   patient-level recognized-mutation recovery.
2. Per-peptide rankings can hide useful mutation-level signal when many routes
   from the same mutation consume the budget.
3. The frozen Epicurus selector is beneficial on Hu_287 and non-regressive on
   Sid, but only two previously inspected patients are evaluable.
4. Epicurus scoring is not generally superior: it beats diversified PRIME by one
   hit on Hu_287 and loses by one hit on Sid.

## What remains unproven

No general superiority claim is justified. Gartner, IMPROVE, multimer, Zhao,
CEDAR, and the Event-B backbone cannot test this exact mechanism from their local
tables because they lack a reliable multi-route mutation-resolved candidate
universe. The next decisive evidence requires reconstructing additional Miller
patients—or another WES/RNA/HLA cohort—then applying this frozen cap-2 crossed
benchmark without changing the policy.
