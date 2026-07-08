import pandas as pd

from epicurus_neo.mhcflurry_features import add_mhcflurry_predictions


class FakePredictor:
    def predict(self, *, peptides, alleles, verbose=0, throw=True):
        return pd.DataFrame(
            {
                "peptide": peptides,
                "affinity": [100.0 + idx for idx, _ in enumerate(peptides)],
                "processing_score": [0.1] * len(peptides),
                "presentation_score": [0.8] * len(peptides),
                "presentation_percentile": [0.5] * len(peptides),
            }
        )


class PairAwareFakePredictor:
    def __init__(self):
        self.calls = []

    def predict(self, *, peptides, alleles, verbose=0, throw=True):
        self.calls.append((tuple(peptides), tuple(alleles)))
        return pd.DataFrame(
            {
                "peptide": peptides,
                "affinity": [10.0] * len(peptides),
                "processing_score": [0.2] * len(peptides),
                "presentation_score": [0.7] * len(peptides),
                "presentation_percentile": [0.4] * len(peptides),
            }
        )


def test_add_mhcflurry_predictions_with_fake_predictor():
    frame = pd.DataFrame(
        {
            "candidate_id": ["a"],
            "source_dataset": ["toy"],
            "study_id": ["s1"],
            "patient_id": ["p1"],
            "hla_allele": ["HLA-A02:01"],
            "mutant_peptide": ["NILGFTFDI"],
            "wildtype_peptide": [""],
            "label": ["positive"],
            "label_weight": [1.0],
            "assay_type": ["synthetic"],
        }
    )
    out = add_mhcflurry_predictions(frame, predictor=FakePredictor())
    assert out.loc[0, "hla_allele_norm"] == "HLA-A*02:01"
    assert out.loc[0, "mhcflurry_affinity"] == 100.0
    assert out.loc[0, "mhcflurry_presentation_score"] == 0.8


def test_add_mhcflurry_predictions_groups_by_allele():
    frame = pd.DataFrame(
        {
            "candidate_id": ["a", "b", "c"],
            "source_dataset": ["toy"] * 3,
            "study_id": ["s1"] * 3,
            "patient_id": ["p1"] * 3,
            "hla_allele": ["HLA-A02:01", "HLA-A02:01", "HLA-B44:02"],
            "mutant_peptide": ["NILGFTFDI", "SIINFEKL", "AEVSVLYTV"],
            "wildtype_peptide": [""] * 3,
            "label": ["positive", "negative", "negative"],
            "label_weight": [1.0] * 3,
            "assay_type": ["synthetic"] * 3,
        }
    )
    predictor = PairAwareFakePredictor()
    out = add_mhcflurry_predictions(frame, predictor=predictor)
    assert len(predictor.calls) == 2
    assert out["mhcflurry_presentation_score"].notna().all()
