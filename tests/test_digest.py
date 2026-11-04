import numpy as np
from src import DigestParams, build_template, simulate
from src.digest import apply_cuts, endonuclease_cuts, exonuclease_trim

# Test Constants
CORE_LENGTH = 147
LINKER_LENGTH = 30


def regular_template(n_chains=50):
    return build_template(50, n_chains, CORE_LENGTH, LINKER_LENGTH, linker_sd=0)


# Test cuts split fragments and ignore cuts outside surviving DNA
def test_apply_cuts():
    starts, ends = np.array([0, 100]), np.array([50, 150])
    cuts = np.array([0, 10, 50, 75, 120, 150])
    starts, ends = apply_cuts(starts, ends, cuts)
    assert list(starts) == [0, 10, 100, 120]
    assert list(ends) == [10, 50, 120, 150]


# Test endonuclease cuts only land in linkers when cores are protected
def test_cuts_avoid_protected_cores():
    template = regular_template()
    cut_rate = np.where(template.seg_is_nuc, 0, 1.0)
    cuts = endonuclease_cuts(template, cut_rate, np.random.default_rng(0))
    is_nuc, _, _ = template.region_at(cuts - 1)
    assert len(cuts) > 0
    assert not np.any(is_nuc)


# Test exonuclease stalls at nucleosome borders
def test_exonuclease_stalls_at_core():
    template = regular_template(n_chains=1)
    params = DigestParams(k_exo_linker=1000, k_exo_nuc=0, dt=1)
    starts, ends = template.chain_start.copy(), template.chain_end.copy()
    starts, ends = exonuclease_trim(template, starts, ends, params, np.random.default_rng(0))
    assert starts[0] == LINKER_LENGTH
    assert ends[0] == template.chain_end[0] - LINKER_LENGTH


# Test total DNA is conserved and fragments stay disjoint without exonuclease
def test_conservation_without_exonuclease():
    template = regular_template()
    params = DigestParams(k_linker=1e-2, k_nuc=1e-4, k_exo_linker=0, k_exo_nuc=0, min_len=0)
    samples = simulate(template, params, [0, 1, 5, 20], np.random.default_rng(0))
    total = template.chain_end[-1]
    for lengths in samples.values():
        assert lengths.sum() == total
    assert len(samples[20]) > len(samples[1]) > len(samples[0]) == 50


# Test fragments stay disjoint and sorted during digestion
def test_fragments_disjoint(monkeypatch):
    import src.digest as digest

    seen = []
    original = digest.exonuclease_trim

    def checked(template, starts, ends, params, rng):
        assert np.all(np.diff(starts) > 0)
        assert np.all(starts[1:] >= ends[:-1])
        seen.append(len(starts))
        return original(template, starts, ends, params, rng)

    monkeypatch.setattr(digest, "exonuclease_trim", checked)
    simulate(regular_template(), DigestParams(), [10], np.random.default_rng(0))
    assert len(seen) > 0


# Test exhaustive digestion of a regular array converges to the core particle
def test_digestion_converges_to_core():
    params = DigestParams(k_linker=1e-2, k_nuc=0, k_exo_linker=5, k_exo_nuc=0)
    samples = simulate(regular_template(), params, [200], np.random.default_rng(0))
    lengths = samples[200]
    assert np.all(lengths % CORE_LENGTH == 0)
    assert np.mean(lengths == CORE_LENGTH) > 0.95


# Test reproducibility with a seed
def test_reproducible():
    template = regular_template()
    a = simulate(template, DigestParams(), [5], np.random.default_rng(1))
    b = simulate(template, DigestParams(), [5], np.random.default_rng(1))
    assert np.array_equal(a[5], b[5])
