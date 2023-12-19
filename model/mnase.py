from model import Fibre
from random import randint, choice
from typing import Literal, Tuple, Optional


class MNase(object):
    @staticmethod
    def choose_fibre(fibres: Tuple[int]) -> int:
        """Chooses a fibre from the given tuple of fibres.

        Args:
            fibres (Tuple[int]): A tuple of integers representing the available fibres.

        Returns:
            int: The chosen fibre.
        """
        return choice(fibres)

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
            Literal[0, 1]: The randomly chosen linker, either 0 (left) or 1 (right).
        """
        return randint(0, 1)

    @staticmethod
    def choose_nth_nucleotide(length: int) -> int:
        """Chooses a random index within the given length.

        Args:
            length (int): The length of the sequence.

        Returns:
            int: The randomly chosen index.
        """
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
            linker (Literal[0, 1]): The chosen linker, either 0 (left) or 1 (right).

        Returns:
            Tuple[Optional[Fibre], Optional[Fibre]]: A tuple containing the
                resulting fibres after cleavage. The first element is the fibre
                before the cleavage position, and the second element is the fibre
                after the cleavage position.

                Special cases:
                    1. nucleosome <= 0      & linker == 0 -> (Fibre, None)
                    2. nucleosome >= length & linker == 1 -> (None, Fibre)
        """
        if (fibre.n_nucleosomes == 1) or _is_terminal_linker(fibre, nucleosome, linker):
            output = _cleave_terminal_linker(fibre, linker)
        else:
            nucleosome = _adjust_nucleosome(fibre, nucleosome, linker)
            output = _cleave_fibre(fibre, nucleosome, linker)
        return output


def _adjust_nucleosome(fibre: Fibre, nucleosome: int, linker: Literal[0, 1]) -> int:
    """Adjust the nucleosome position.

    Args:
        fibre (Fibre): The fibre to be cleaved.
        nucleosome (int): The position of the nucleosome.
        linker (Literal[0, 1]): The chosen linker, either 0 (left) or 1 (right).

    Returns:
        int: The adjusted nucleosome position.
    """
    if nucleosome <= 0 and linker == 1:
        adjusted_nucleosome = nucleosome + 1
    elif nucleosome >= fibre.n_nucleosomes and linker == 0:
        adjusted_nucleosome = nucleosome - 1
    else:
        adjusted_nucleosome = nucleosome
    return adjusted_nucleosome


def _is_terminal_linker(fibre: Fibre, nucleosome: int, linker: Literal[0, 1]) -> bool:
    """Determines if a linker is a terminal linker in a given nucleosome.

    Args:
        fibre (Fibre): The fibre object.
        nucleosome (int): The index of the nucleosome.
        linker (Literal[0, 1]): The chosen linker, either 0 (left) or 1 (right).

    Returns:
        bool: True if the linker is a terminal linker, False otherwise.
    """
    return _is_terminal_left_linker(nucleosome, linker) or _is_terminal_right_linker(
        fibre, nucleosome, linker
    )


def _is_terminal_left_linker(nucleosome: int, linker: Literal[0, 1]) -> bool:
    """Checks if the given nucleosome is a terminal left linker.

    Args:
        nucleosome (int): The index of the nucleosome.
        linker (Literal[0, 1]): The chosen linker, either 0 (left) or 1 (right).

    Returns:
        bool: True if the nucleosome is a terminal left linker, False otherwise.
    """
    return (nucleosome <= 0) and (linker == 0)


def _is_terminal_right_linker(
    fibre: Fibre, nucleosome: int, linker: Literal[0, 1]
) -> bool:
    """Checks if the given nucleosome is the last one in the fibre and if the linker is on the right side.

    Args:
        fibre (Fibre): The fibre object.
        nucleosome (int): The index of the nucleosome.
        linker (Literal[0, 1]): The chosen linker, either 0 (left) or 1 (right).

    Returns:
        bool: True if the nucleosome is the last one in the fibre and the linker
            is on the right side, False otherwise.
    """
    return (nucleosome >= fibre.n_nucleosomes) and (linker == 1)


def _cleave_fibre(
    fibre: Fibre, nucleosomes: int, linker: Literal[0, 1]
) -> Tuple[Fibre, Fibre]:
    """Cleave a fibre into two fibres based on the given parameters.

    Args:
        fibre (Fibre): The fibre to be cleaved.
        nucleosomes (int): The number of nucleosomes.
        linker (Literal[0, 1]): The chosen linker, either 0 (left) or 1 (right).

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
        linker (Literal[0, 1]): The chosen linker, either 0 (left) or 1 (right).

    Returns:
        Tuple[Optional[Fibre], Optional[Fibre]]: A tuple containing the modified fibre(s).
        If linker is 0, the left linker is cleaved and the modified fibre is returned as the first element of the tuple.
        If linker is 1, the right linker is cleaved and the modified fibre is returned as the second element of the tuple.
    """
    output = [None, None]
    linker_length = fibre.get_linker(linker)
    if linker_length <= 0:
        output[linker] = fibre
    else:
        n_nucleotides_removed = MNase.remove_n_nucleotides(linker_length)
        if linker_length < n_nucleotides_removed:
            n_nucleotides_removed = linker_length

        fibre.set_linker(linker, linker_length - n_nucleotides_removed)
        output[linker] = fibre
    return tuple(output)
