from model import Fibre
from random import randint
from typing import Literal, Tuple, Optional, List


class MNase(object):
    @staticmethod
    def choose_fibre(fibres: List[Fibre]) -> int:
        """
        Chooses a random fibre from the given list of fibres.

        Args:
            fibres (List[Fibre]): A list of fibres to choose from.

        Returns:
            int: The index of the chosen fibre.

        Raises:
            IndexError: If the list of fibres is empty.
        """
        remove_nones = [fibre for fibre in fibres if fibre is not None]
        return randint(0, len(remove_nones) - 1)

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
    def choose_nth_nucleotide(length: int) -> int:
        return randint(0, length)

    @staticmethod
    def remove_n_nucleotides(length: int) -> int:
        """Removes a random number of nucleotides from a given length.

        Args:
            length (int): The length of the nucleotide sequence.

        Returns:
            int: The number of nucleotides to be removed.
        """
        return randint(1, length)

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
                    1. nucleosome <= 0      & linker == 1 -> (None, Fibre)
                    2. nucleosome >= length & linker == 0 -> (Fibre, None)
        """
        if nucleosome <= 0 and linker == 1:
            output = _cleave_terminal_linker(fibre, linker)
        elif nucleosome >= fibre.n_nucleosomes and linker == 0:
            output = _cleave_terminal_linker(fibre, linker)
        elif fibre.n_nucleosomes == 1:
            output = _cleave_terminal_linker(fibre, linker)
        else:
            output = _cleave_fibre(fibre, nucleosome, linker)

        return output


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
    class_constants = {
        "linker_length": fibre.LINKER_LENGTH,
        "nucleosome_length": fibre.NUCLEOSOME_LENGTH,
    }

    cleavage_position = MNase.choose_nth_nucleotide(fibre.get_linker(linker))

    fibre_before = Fibre(nucleosomes, **class_constants)
    fibre_before.set_linker(linker, cleavage_position)

    fibre_after = Fibre(fibre.n_nucleosomes - nucleosomes, **class_constants)
    fibre_after.set_opposing_linker(
        linker, fibre.get_linker(linker) - cleavage_position
    )

    return (fibre_before, fibre_after)


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
