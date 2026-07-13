# Frozen portfolio generalization stress test

> Two previously inspected patients; descriptive and post-hoc. Hu_287 is the discovery replay, Sid is a stress test. This is not independent validation.

## Primary result — k=20, frozen cap=2

| Patient | PRIME plain | PRIME + selector | Epicurus plain | Epicurus + selector | Selector delta on PRIME | Selector delta on Epicurus |
|---|---:|---:|---:|---:|---:|---:|
| Hu_287 | 0/3 | 2/3 | 0/3 | 3/3 | +2 | +3 |
| Sid | 2/3 | 2/3 | 1/3 | 1/3 | +0 | +0 |

## Slot-use diagnostics

| Patient | Arm | Slots | Unique mutations | Duplicate burden |
|---|---|---:|---:|---:|
| Hu_287 | `prime_plain` | 20 | 8 | 12 |
| Hu_287 | `prime_route_aware` | 20 | 11 | 9 |
| Hu_287 | `epicurus_plain` | 20 | 3 | 17 |
| Hu_287 | `epicurus_route_aware` | 20 | 10 | 10 |
| Sid | `prime_plain` | 20 | 17 | 3 |
| Sid | `prime_route_aware` | 20 | 17 | 3 |
| Sid | `epicurus_plain` | 20 | 8 | 12 |
| Sid | `epicurus_route_aware` | 20 | 13 | 7 |

## Interpretation guardrail

The crossed control is decisive for attribution: if PRIME gains similarly from the selector, the supported mechanism is scorer-agnostic portfolio diversification. Sid determines whether the Hu_287 mechanism survives a second mutation-resolved patient; a tie or regression must be reported as such.

Sensitivity results and the full eligibility audit are in `RESULT.json`. Caps other than 2 and k values other than 20 are post-hoc diagnostics, not alternative headline endpoints.
