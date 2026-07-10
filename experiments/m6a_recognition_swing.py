"""Thin CLI entry point for the M6A Event-B-only recognition swing.

The reusable implementation lives in :mod:`epicurus_neo.m6.runner`.
Run with: ``.venv/bin/python experiments/m6a_recognition_swing.py``
"""

from __future__ import annotations

from epicurus_neo.m6.runner import run

if __name__ == "__main__":
    result = run()
    print(f"universal verdict: {result['universal']['verdict']}")
    print(f"presentation verdict: {result['presentation']['verdict']}")
