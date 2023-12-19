from typing import Literal
from model.validators import validate_fibre_input


class Fibre(object):
    def __init__(self, n_nucleosomes: int, linker_length: int, nucleosome_length: int):
        validate_fibre_input(n_nucleosomes, linker_length, nucleosome_length)

        # Constants
        self.LINKER_LENGTH = linker_length
        self.NUCLEOSOME_LENGTH = nucleosome_length

        # Attributes
        self.n_nucleosomes = n_nucleosomes
        self.left_linker = linker_length
        self.right_linker = linker_length

    def n_nucleotides(self) -> int:
        """Calculates the total number of nucleotides in the fibre.

        Args:
            nucleosome_dna_length (int): The length of the DNA in a nucleosome.

        Returns:
            int: The total number of nucleotides in the fibre.
        """
        nucleosome_length = self.NUCLEOSOME_LENGTH * self.n_nucleosomes
        terminal_linker_length = self.left_linker + self.right_linker
        inner_linker_length = self.LINKER_LENGTH * (self.n_nucleosomes - 1)
        return nucleosome_length + terminal_linker_length + inner_linker_length

    def get_linker(self, index: Literal[0, 1]) -> int:
        """Returns the linker value at the specified index.

        Args:
            index (Literal[0, 1]): The linker index, either 0 (left) or 1 (right).

        Returns:
            int: The linker value at the specified index.

        Raises:
            ValueError: If the index is not 0 or 1.
        """
        if index == 0:
            return self.left_linker
        elif index == 1:
            return self.right_linker
        else:
            raise ValueError("Invalid index. Must be either 0 or 1.")

    def set_linker(self, index: Literal[0, 1], length: int) -> None:
        """Sets the length of the linker for the specified index.

        Args:
            index (Literal[0, 1]): The linker index, either 0 (left) or 1 (right).
            length (int): The length of the linker to set.

        Raises:
            ValueError: If the index is not 0 or 1.
        """
        if index == 0:
            self.left_linker = length
        elif index == 1:
            self.right_linker = length
        else:
            raise ValueError("Invalid index. Must be either 0 or 1.")

    def set_opposing_linker(self, index: Literal[0, 1], length: int) -> None:
        """
        Sets the length of the opposing linker based on the given index.

        Args:
            index (Literal[0, 1]): The linker index, either 0 (left) or 1 (right).
            length (int): The length of the linker.

        Raises:
            ValueError: If the index is not 0 or 1.
        """
        if index == 0:
            self.right_linker = length
        elif index == 1:
            self.left_linker = length
        else:
            raise ValueError("Invalid index. Must be either 0 or 1.")

    def __repr__(self):
        return f"ID: {id(self)}\tLeft: {self.left_linker}\tRight: {self.right_linker}\tNucleosomes: {self.n_nucleosomes}"
