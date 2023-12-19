import pytest
from model import Fibre

# Test Constants
VALID_N_NUCLEOSOMES = 10
VALID_LINKER_LENGTH = 10
VALID_NUCLEOSOME_LENGTH = 147


# Test instantiation
def test_fibre_instantiation():
    fibre = Fibre(VALID_N_NUCLEOSOMES, VALID_LINKER_LENGTH, VALID_NUCLEOSOME_LENGTH)
    assert fibre.n_nucleosomes == VALID_N_NUCLEOSOMES
    assert fibre.LINKER_LENGTH == VALID_LINKER_LENGTH
    assert fibre.NUCLEOSOME_LENGTH == VALID_NUCLEOSOME_LENGTH


# Test invalid instantiation
@pytest.mark.parametrize(
    "n_nucleosomes, linker_length, nucleosome_length",
    [
        (-1, VALID_LINKER_LENGTH, VALID_NUCLEOSOME_LENGTH),  # Invalid nucleosomes
        (VALID_N_NUCLEOSOMES, -1, VALID_NUCLEOSOME_LENGTH),  # Invalid linker length
        (VALID_N_NUCLEOSOMES, VALID_LINKER_LENGTH, -1),  # Invalid nucleosome length
    ],
)
def test_fibre_invalid_instantiation(n_nucleosomes, linker_length, nucleosome_length):
    with pytest.raises(ValueError):
        Fibre(n_nucleosomes, linker_length, nucleosome_length)


# Test for invalid type for linker length
def test_invalid_type_linker_length():
    with pytest.raises(TypeError):
        Fibre(VALID_N_NUCLEOSOMES, "invalid_type", VALID_NUCLEOSOME_LENGTH)


# Test n_nucleotides calculation
def test_n_nucleotides():
    fibre = Fibre(VALID_N_NUCLEOSOMES, VALID_LINKER_LENGTH, VALID_NUCLEOSOME_LENGTH)
    expected_length = (
        VALID_NUCLEOSOME_LENGTH * VALID_N_NUCLEOSOMES
        + 2 * VALID_LINKER_LENGTH
        + VALID_LINKER_LENGTH * (VALID_N_NUCLEOSOMES - 1)
    )
    assert fibre.n_nucleotides() == expected_length


# Test get_linker method
@pytest.mark.parametrize(
    "index, expected_length",
    [(0, VALID_LINKER_LENGTH), (1, VALID_LINKER_LENGTH)],  # Left linker  # Right linker
)
def test_get_linker(index, expected_length):
    fibre = Fibre(VALID_N_NUCLEOSOMES, VALID_LINKER_LENGTH, VALID_NUCLEOSOME_LENGTH)
    assert fibre.get_linker(index) == expected_length


# Test get_linker with invalid index
def test_get_linker_invalid_index():
    fibre = Fibre(VALID_N_NUCLEOSOMES, VALID_LINKER_LENGTH, VALID_NUCLEOSOME_LENGTH)
    with pytest.raises(ValueError):
        fibre.get_linker(2)  # Invalid index


# Test set_linker method
def test_set_linker():
    fibre = Fibre(VALID_N_NUCLEOSOMES, VALID_LINKER_LENGTH, VALID_NUCLEOSOME_LENGTH)
    new_length = 20
    fibre.set_linker(0, new_length)
    assert fibre.left_linker == new_length
    fibre.set_linker(1, new_length)
    assert fibre.right_linker == new_length


# Test set_linker with invalid index
def test_set_linker_invalid_index():
    fibre = Fibre(VALID_N_NUCLEOSOMES, VALID_LINKER_LENGTH, VALID_NUCLEOSOME_LENGTH)
    with pytest.raises(ValueError):
        fibre.set_linker(2, 20)  # Invalid index


# Test set_opposing_linker method
def test_set_opposing_linker():
    fibre = Fibre(VALID_N_NUCLEOSOMES, VALID_LINKER_LENGTH, VALID_NUCLEOSOME_LENGTH)
    new_length = 20
    fibre.set_opposing_linker(0, new_length)
    assert fibre.right_linker == new_length
    fibre.set_opposing_linker(1, new_length)
    assert fibre.left_linker == new_length


# Test set_opposing_linker with invalid index
def test_set_opposing_linker_invalid_index():
    fibre = Fibre(VALID_N_NUCLEOSOMES, VALID_LINKER_LENGTH, VALID_NUCLEOSOME_LENGTH)
    with pytest.raises(ValueError):
        fibre.set_opposing_linker(2, 20)  # Invalid index


# Test for zero-length linkers
def test_zero_length_linkers():
    fibre = Fibre(VALID_N_NUCLEOSOMES, 0, VALID_NUCLEOSOME_LENGTH)
    assert fibre.n_nucleotides() == VALID_NUCLEOSOME_LENGTH * VALID_N_NUCLEOSOMES
