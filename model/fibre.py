class Fibre(object):
    def __init__(self, n_nucleosomes: int, linker_length: int):
        self.n_nucleosomes = n_nucleosomes
        self.left_linker = linker_length
        self.right_linker = linker_length

    def n_nucleotides(self, nucleosome_dna_length: int) -> int:
        """Calculates the total number of nucleotides in the fibre.

        Args:
            nucleosome_dna_length (int): The length of the DNA in a nucleosome.

        Returns:
            int: The total number of nucleotides in the fibre.
        """
        total_nucleosome_dna_length = nucleosome_dna_length * self.n_nucleosomes
        linker_dna_length = self.left_linker + self.right_linker
        return total_nucleosome_dna_length + linker_dna_length
