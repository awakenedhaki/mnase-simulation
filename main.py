from tqdm import tqdm
from model import Fibre, MNase


def main():
    TIME_STEPS = 10_000
    N_NUCLEOSOMES = 10_000_000
    NUCLEOSOME_LENGTH = 147
    LINKER_LENGTH = 10

    # Initialize fibre array
    fibre = Fibre(N_NUCLEOSOMES, LINKER_LENGTH, NUCLEOSOME_LENGTH)
    fibres = {id(fibre): fibre}

    # Initialize model
    for _ in tqdm(range(time_steps)):
        # Select which fibre to cleave
        selected_fibre = MNase.choose_fibre(tuple(fibres.keys()))
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
        if all(fibre is not None for fibre in cleaved_products):
            fibre_before = cleaved_products[0]
            fibre_after = cleaved_products[1]
            fibres[id(fibre_before)] = fibre_before
            fibres[id(fibre_after)] = fibre_after

            if not id(fibre_before) == selected_fibre:
                del fibres[id(fibre_before)]
            elif not id(fibre_after) == selected_fibre:
                del fibres[id(fibre_after)]
        elif fibre := next(filter(lambda x: x is not None, cleaved_products)):
            fibres[selected_fibre] = fibre

    fibre_lengths = [fibre.n_nucleotides() for fibre in fibres.values()]
    return fibre_lengths


if __name__ == "__main__":
    print("\n", main())
