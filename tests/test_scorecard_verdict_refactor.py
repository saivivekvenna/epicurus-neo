from benchmark.scorecard import pre_registered_verdict


def test_verdict_rule_matches_registered_semantics():
    accept = {"delta_vs_baseline": 0.5, "delta_ci": [0.1, 0.9]}
    safe = {"delta_vs_baseline": 0.0, "delta_ci": [-0.1, 0.2]}
    assert pre_registered_verdict(accept, safe, safe) == "ACCEPT"
    unresolved = {"delta_vs_baseline": 0.3, "delta_ci": [-0.2, 0.8]}
    assert pre_registered_verdict(unresolved, safe, safe) == "CONSISTENT_WITH_NO_EFFECT"
    regressed = {"delta_vs_baseline": -0.4, "delta_ci": [-0.9, -0.1]}
    assert pre_registered_verdict(regressed, safe, safe) == "REJECT"
    # A significant primary gain is vetoed by a significant co-primary regression.
    harm = {"delta_vs_baseline": -0.3, "delta_ci": [-0.5, -0.1]}
    assert pre_registered_verdict(accept, harm, safe) == "REJECT"
