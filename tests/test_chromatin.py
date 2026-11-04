import numpy as np
import pytest
from src import build_template

# Test Constants
VALID_N_NUCLEOSOMES = 10
VALID_N_CHAINS = 3
VALID_CORE_LENGTH = 147
VALID_LINKER_LENGTH = 20


def regular_template():
    return build_template(
        VALID_N_NUCLEOSOMES,
        VALID_N_CHAINS,
        VALID_CORE_LENGTH,
        linker_mean=VALID_LINKER_LENGTH,
        linker_sd=0,
    )


# Test segments are contiguous and cover every chain
def test_segments_contiguous():
    template = regular_template()
    assert np.all(template.seg_start[1:] == template.seg_end[:-1])
    assert template.seg_start[0] == 0
    assert template.seg_end[-1] == template.chain_end[-1]
    assert np.all(template.chain_start[1:] == template.chain_end[:-1])


# Test regular array lengths
def test_regular_template_lengths():
    template = regular_template()
    expected_chain = (
        VALID_CORE_LENGTH * VALID_N_NUCLEOSOMES
        + VALID_LINKER_LENGTH * (VALID_N_NUCLEOSOMES + 1)
    )
    assert np.all(template.chain_end - template.chain_start == expected_chain)
    assert template.seg_is_nuc.sum() == VALID_N_NUCLEOSOMES * VALID_N_CHAINS
    assert np.all(template.seg_length[template.seg_is_nuc] == VALID_CORE_LENGTH)
    assert np.all(template.seg_length[~template.seg_is_nuc] == VALID_LINKER_LENGTH)


# Test variable linkers are reproducible and non-negative
def test_variable_linkers():
    kwargs = dict(linker_mean=VALID_LINKER_LENGTH, linker_sd=15)
    a = build_template(1000, rng=np.random.default_rng(0), **kwargs)
    b = build_template(1000, rng=np.random.default_rng(0), **kwargs)
    linkers = a.seg_length[~a.seg_is_nuc]
    assert np.array_equal(a.seg_end, b.seg_end)
    assert np.all(linkers > 0)
    assert linkers.std() > 0


# Test zero-length linkers are dropped
def test_zero_length_linkers():
    template = build_template(VALID_N_NUCLEOSOMES, linker_mean=0)
    assert np.all(template.seg_is_nuc)
    assert template.chain_end[0] == VALID_CORE_LENGTH * VALID_N_NUCLEOSOMES


# Test nucleosome-free regions
def test_nucleosome_free_regions():
    template = build_template(VALID_N_NUCLEOSOMES, nfr_prob=1, nfr_length=500)
    assert np.all(template.seg_length[~template.seg_is_nuc] == 500)


# Test region lookup
def test_region_at():
    template = regular_template()
    is_nuc, start, end = template.region_at(np.array([0, VALID_LINKER_LENGTH]))
    assert list(is_nuc) == [False, True]
    assert list(start) == [0, VALID_LINKER_LENGTH]
    assert list(end) == [VALID_LINKER_LENGTH, VALID_LINKER_LENGTH + VALID_CORE_LENGTH]


# Test invalid instantiation
@pytest.mark.parametrize(
    "kwargs",
    [
        dict(n_nucleosomes=0),
        dict(n_nucleosomes=VALID_N_NUCLEOSOMES, n_chains=0),
        dict(n_nucleosomes=VALID_N_NUCLEOSOMES, linker_mean=-1),
        dict(n_nucleosomes=VALID_N_NUCLEOSOMES, linker_sd=-1),
        dict(n_nucleosomes=VALID_N_NUCLEOSOMES, nfr_prob=1.5),
    ],
)
def test_invalid_input(kwargs):
    with pytest.raises(ValueError):
        build_template(**kwargs)


# Test for invalid type
def test_invalid_type():
    with pytest.raises(TypeError):
        build_template("invalid_type")
