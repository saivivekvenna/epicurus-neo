# Auto-Research Workflow

Auto-research is used as an experiment generator and error analyst. It is not a
direct peptide-scoring model.

## Loop

1. Normalize public data.
2. Run leakage-aware grouped CV or a frozen train/test evaluation.
3. Save scored candidates.
4. Generate a failure report and LLM hypothesis prompt.
5. Ask an LLM to propose testable feature/model changes in the required YAML
   format.
6. Implement one hypothesis behind a config or feature flag.
7. Re-run the benchmark.
8. Keep the hypothesis only if it improves held-out top-k metrics.

## Generate Research Artifacts

```bash
epicurus research-report \
  --scored outputs/validation.scored.csv \
  --group-col patient_id \
  --score-col epicurus_score \
  -k 20 \
  --output-dir outputs/research/run_001
```

Outputs:

```text
outputs/research/run_001/failure_report.json
outputs/research/run_001/hypothesis_prompt.md
```

## Rules

- The LLM may propose hypotheses.
- The benchmark accepts or rejects hypotheses.
- No LLM direct peptide scoring.
- No locked test labels in hypothesis search.
- Improvements must be measured on grouped, leakage-controlled validation.

