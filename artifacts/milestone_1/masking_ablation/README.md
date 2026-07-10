# Ten-seed masking ablation question sets

Generated from the official IMPROVE CV table in the cloned `SRHgroup/IMPROVE_paper` repository by:

```bash
python scripts/milestone_1.py /path/to/IMPROVE_paper generate-ablation \
  artifacts/milestone_1/masking_ablation
```

Each seed has 50 balanced questions per condition (25 positive, 25 negative), with A and B disjoint
within that seed. `questions/` contains no labels. `answer_keys/` must remain hidden until all LLM
answers for that seed are committed. Condition A exposes gene, protein-position/change, mutation
consequence, and HLA. Condition B exposes mutant peptide, wild-type peptide, and HLA, with the gene
and mutation identifier masked.

No LLM batch inference endpoint is available in this workspace, so this milestone emits the ten
question sets as required and does not fabricate replication results.
