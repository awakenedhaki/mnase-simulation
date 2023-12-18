from model import Fibre
from random import randint
from typing import Literal, Tuple, Optional


class MNase(object):
    @staticmethod
    def choose_nucleosome(length: int) -> int:
        """Randomly chooses a nucleosome index within the given length.

        Args:
            length (int): The length of the fibre.

        Returns:
            int: The randomly chosen nucleosome index.
        """
        return randint(0, length)

    @staticmethod
    def choose_linker() -> Literal[0, 1]:
        """Randomly chooses a linker.

        Returns:
            Literal[0, 1]: The randomly chosen linker, either 0 (right) or 1 (left).
        """
        return randint(0, 1)

    @staticmethod
    def cleave_fibre(
        fibre: Fibre, nucleosome: int, linker: Literal[0, 1]
    ) -> Tuple[Optional[Fibre], Optional[Fibre]]:
        """Cleaves a fibre at a specific nucleosome position.

        Args:
            fibre (Fibre): The fibre to be cleaved.
            nucleosome (int): The position of the nucleosome to cleave at.
            linker (Literal[0, 1]): The chosen linker, either 0 (right) or 1 (left).

        Returns:
            Tuple[Optional[Fibre], Optional[Fibre]]: A tuple containing the
                resulting fibres after cleavage. The first element is the fibre
                before the cleavage position, and the second element is the fibre
                after the cleavage position.

                Special cases:
                    1. nucleosome == 0      & linker == 0 -> (None, Fibre)
                    2. nucleosome == length & linker == 1 -> (Fibre, None)
        """
        if nucleosome <= 0 and linker == 0:
            output = _cleave_terminal_linker(fibre, linker)
        elif nucleosome >= fibre.n_nucleosomes and linker == 1:
            output = _cleave_terminal_linker(fibre, linker)
        else:
            output = _cleave_fibre(fibre, nucleosome, linker)

        return output


def _remove_n_nucleotides(length: int) -> int:
    """Removes a random number of nucleotides from a given length.

    Args:
        length (int): The length of the nucleotide sequence.

    Returns:
        int: The number of nucleotides to be removed.
    """
    return randint(1, length)


def _cleave_fibre(
    fibre: Fibre, nucleosomes: int, linker: Literal[0, 1]
) -> Tuple[Fibre, Fibre]:
    """Cleave a fibre into two fibres based on the given parameters.

    Args:
        fibre (Fibre): The fibre to be cleaved.
        nucleosomes (int): The number of nucleosomes.
        linker (Literal[0, 1]): The chosen linker, either 0 (right) or 1 (left).

    Returns:
        Tuple[Fibre, Fibre]: A tuple containing the two resulting fibres.
    """
    linker_length, n_nucleosomes = fibre.LINKER_LENGTH, fibre.n_nucleosomes

    # Calculating linker length
    linker_length_1 = randint(1, linker_length)
    linker_length_2 = linker_length - linker_length_1

    shortest_length = min(linker_length_1, linker_length_2)
    longest_length = max(linker_length_1, linker_length_2)

    # Calculating nucleosome number
    nucleosomes_1 = nucleosomes
    nucleosomes_2 = n_nucleosomes - nucleosomes

    lowest_nucleosomes = min(nucleosomes_1, nucleosomes_2)
    highest_nucleosomes = max(nucleosomes_1, nucleosomes_2)

    # Make output fibres
    fibres = [None, None]
    if linker == 0:
        fibres[0] = Fibre(lowest_nucleosomes, linker_length, fibre.NUCLEOSOME_LENGTH)
        fibres[0].right_linker = shortest_length
        fibres[1] = Fibre(highest_nucleosomes, linker_length, fibre.NUCLEOSOME_LENGTH)
        fibres[1].left_linker = longest_length
    elif linker == 1:
        fibres[0] = Fibre(highest_nucleosomes, linker_length, fibre.NUCLEOSOME_LENGTH)
        fibres[0].right_linker = longest_length
        fibres[1] = Fibre(lowest_nucleosomes, linker_length, fibre.NUCLEOSOME_LENGTH)
        fibres[1].left_linker = shortest_length

    return tuple(fibres)


def _cleave_terminal_linker(
    fibre: Fibre, linker: Literal[0, 1]
) -> Tuple[Optional[Fibre], Optional[Fibre]]:
    """Cleave the terminal linker of a fibre.

    Args:
        fibre (Fibre): The fibre to cleave the linker from.
        linker (Literal[0, 1]): The chosen linker, either 0 (right) or 1 (left).

    Returns:
        Tuple[Optional[Fibre], Optional[Fibre]]: A tuple containing the modified fibre(s).
        If linker is 0, the right linker is cleaved and the modified fibre is returned as the second element of the tuple.
        If linker is 1, the left linker is cleaved and the modified fibre is returned as the first element of the tuple.
    """
    linker_length = fibre.left_linker if linker == 0 else fibre.right_linker
    n_nucleotides_removed = _remove_n_nucleotides(linker_length)

    output = (None, None)
    if linker == 0:
        fibre.right_linker -= n_nucleotides_removed
        output = (None, fibre)
    else:
        fibre.left_linker -= n_nucleotides_removed
        output = (fibre, None)

    return output
