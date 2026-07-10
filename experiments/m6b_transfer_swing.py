"""Thin CLI entry point for the M6B Event-A -> Event-B transfer swing.

The reusable implementation lives in :mod:`epicurus_neo.m6.transfer_runner`.
Run with: ``.venv/bin/python experiments/m6b_transfer_swing.py``
"""

from __future__ import annotations

from epicurus_neo.m6.transfer_runner import run_m6b

if __name__ == "__main__":
    result = run_m6b()
    transfer = result["transfer"]
    print(f"transfer verdict: {transfer['verdict']}")
    print(f"macro AUROC delta: {transfer['macro_auroc_delta']:.4f}")
    print(f"folds improved: {transfer['folds_improved']}/{transfer['n_folds_scored']}")
