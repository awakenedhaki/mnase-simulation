from model import Fibre, MNase


def main():
    TIME_STEPS = 25
    N_NUCLEOSOMES = 1000
    NUCLEOSOME_LENGTH = 147
    LINKER_LENGTH = 10

    # Initialize fibre array
    fibres = [None] * N_NUCLEOSOMES
    fibres[0] = Fibre(N_NUCLEOSOMES, LINKER_LENGTH, NUCLEOSOME_LENGTH)

    # Initialize model
    for i, _ in enumerate(range(TIME_STEPS)):
        # Select which fibre to cleave
        selected_fibre = MNase.choose_fibre(fibres)
        fibre = fibres[selected_fibre]

        # Select which nucleosome to cleave
        selected_nucleosome = MNase.choose_nucleosome(fibre.n_nucleosomes)

        # Select which linker to cleave
        selected_linker = MNase.choose_linker()

        # Cleave fibre
        cleaved_products = MNase.cleave_fibre(
            fibre, selected_nucleosome, selected_linker
        )

        # Update fibre array
        if any(fibre is not None for fibre in cleaved_products) and any(
            fibre is None for fibre in fibres
        ):
            fibres[selected_fibre] = cleaved_products[0]
            next_none_element_idx = next(
                i for i, fibre in enumerate(fibres) if fibre is None
            )
            fibres[next_none_element_idx] = cleaved_products[1]
        elif fibre := next(filter(lambda x: x is not None, cleaved_products)):
            fibres[selected_fibre] = fibre

    fibre_lengths = [fibre.n_nucleotides() for fibre in fibres if fibre is not None]
    return fibre_lengths


if __name__ == "__main__":
    print("\n", main())
