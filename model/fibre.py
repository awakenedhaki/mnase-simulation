class Fibre(object):
    def __init__(self, n_nucleosomes: int, linker_length: int, nucleosome_length: int):
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
        total_nucleosome_dna_length = self.NUCLEOSOME_LENGTH * self.n_nucleosomes
        linker_dna_length = self.left_linker + self.right_linker
        return total_nucleosome_dna_length + linker_dna_length

    def __repr__(self):
        return f"Fibre({self.n_nucleosomes}, {self.LINKER_LENGTH}, {self.NUCLEOSOME_LENGTH})"
