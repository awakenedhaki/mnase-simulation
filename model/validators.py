def validate_fibre_input(
    n_nucleosomes: int, linker_length: int, nucleosome_length: int
):
    """Validate inputs for the Fibre class.

    Checks if the inputs are integers and positive. Raises appropriate errors
    if the validation fails.

    Args:
        n_nucleosomes: An integer representing the number of nucleosomes.
        linker_length: An integer representing the length of the linker DNA.
        nucleosome_length: An integer representing the length of the nucleosome DNA.

    Raises:
        TypeError: If any of the inputs is not an integer.
        ValueError: If any of the inputs is not positive.
    """
    if not all(
        isinstance(value, int)
        for value in [n_nucleosomes, linker_length, nucleosome_length]
    ):
        raise TypeError(
            "n_nucleosomes, linker_length, and nucleosome_length must be integers"
        )

    if not all(
        value >= 0 for value in [n_nucleosomes, linker_length, nucleosome_length]
    ):
        raise ValueError(
            "n_nucleosomes, linker_length, and nucleosome_length must be positive"
        )
