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
            ...
        elif nucleosome >= fibre.n_nucleosomes and linker == 1:
            ...
        else:
            ...
