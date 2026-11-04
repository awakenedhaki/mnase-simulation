from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class Template:
    """Linear chromatin map made of alternating linker and nucleosome segments.

    Segments are half-open intervals [seg_start, seg_end) in bp coordinates,
    sorted and contiguous. Several independent chains (DNA molecules) are laid
    end to end; `chain_start` / `chain_end` mark the boundaries of each molecule.

    Attributes:
        seg_start (np.ndarray): Start coordinate of every segment.
        seg_end (np.ndarray): End coordinate (exclusive) of every segment.
        seg_is_nuc (np.ndarray): True where the segment is a nucleosome core.
        chain_start (np.ndarray): Start coordinate of every chain.
        chain_end (np.ndarray): End coordinate (exclusive) of every chain.
    """

    seg_start: np.ndarray
    seg_end: np.ndarray
    seg_is_nuc: np.ndarray
    chain_start: np.ndarray
    chain_end: np.ndarray

    @property
    def seg_length(self) -> np.ndarray:
        return self.seg_end - self.seg_start

    def segment_at(self, pos: np.ndarray) -> np.ndarray:
        """Returns the index of the segment containing each bp position.

        Args:
            pos (np.ndarray): bp coordinates (must lie inside the template).

        Returns:
            np.ndarray: Segment index for every position.
        """
        return np.searchsorted(self.seg_start, pos, side="right") - 1

    def region_at(self, pos: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Returns whether each position is nucleosomal, plus its segment bounds.

        Args:
            pos (np.ndarray): bp coordinates (must lie inside the template).

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]: (is_nuc, seg_start, seg_end)
                for the segment containing each position.
        """
        idx = self.segment_at(pos)
        return self.seg_is_nuc[idx], self.seg_start[idx], self.seg_end[idx]


def validate_template_input(
    n_nucleosomes: int,
    n_chains: int,
    core_length: int,
    linker_mean: float,
    linker_sd: float,
    nfr_prob: float,
    nfr_length: int,
) -> None:
    """Validate inputs for build_template.

    Raises:
        TypeError: If a count or length is not an integer.
        ValueError: If a value is out of range.
    """
    integers = [n_nucleosomes, n_chains, core_length, nfr_length]
    if not all(isinstance(value, (int, np.integer)) for value in integers):
        raise TypeError(
            "n_nucleosomes, n_chains, core_length and nfr_length must be integers"
        )
    if n_nucleosomes < 1 or n_chains < 1 or core_length < 1:
        raise ValueError("n_nucleosomes, n_chains and core_length must be >= 1")
    if linker_mean < 0 or linker_sd < 0 or nfr_length < 0:
        raise ValueError("linker_mean, linker_sd and nfr_length must be >= 0")
    if not 0 <= nfr_prob <= 1:
        raise ValueError("nfr_prob must be between 0 and 1")


def build_template(
    n_nucleosomes: int,
    n_chains: int = 1,
    core_length: int = 147,
    linker_mean: float = 30,
    linker_sd: float = 0,
    nfr_prob: float = 0,
    nfr_length: int = 300,
    rng: Optional[np.random.Generator] = None,
) -> Template:
    """Builds a chromatin template of independent nucleosome arrays.

    Each chain is laid out as linker, core, linker, core, ..., core, linker.
    Linker lengths are drawn from a normal distribution (clipped at 0 and
    rounded); with linker_sd=0 the array is perfectly regular. Each linker is
    replaced by a nucleosome-free region of nfr_length bp with probability nfr_prob.

    Args:
        n_nucleosomes (int): Nucleosomes per chain.
        n_chains (int): Number of independent DNA molecules.
        core_length (int): DNA wrapped around each nucleosome core (bp).
        linker_mean (float): Mean linker length (bp).
        linker_sd (float): Standard deviation of linker length (bp).
        nfr_prob (float): Probability that a linker is a nucleosome-free region.
        nfr_length (int): Length of a nucleosome-free region (bp).
        rng (Optional[np.random.Generator]): Random number generator.

    Returns:
        Template: The chromatin map.
    """
    validate_template_input(
        n_nucleosomes, n_chains, core_length, linker_mean, linker_sd, nfr_prob, nfr_length
    )
    rng = np.random.default_rng() if rng is None else rng

    n_linkers = n_nucleosomes + 1
    linkers = rng.normal(linker_mean, linker_sd, size=(n_chains, n_linkers))
    linkers = np.clip(np.rint(linkers), 0, None).astype(np.int64)
    if nfr_prob > 0:
        linkers[rng.random(linkers.shape) < nfr_prob] = nfr_length

    # Interleave linkers (even columns) with cores (odd columns)
    lengths = np.empty((n_chains, 2 * n_nucleosomes + 1), dtype=np.int64)
    lengths[:, 0::2] = linkers
    lengths[:, 1::2] = core_length
    is_nuc = np.zeros(lengths.shape, dtype=bool)
    is_nuc[:, 1::2] = True

    # Drop zero-length linkers so every segment covers at least one bp
    lengths, is_nuc = lengths.ravel(), is_nuc.ravel()
    chain_length = lengths.reshape(n_chains, -1).sum(axis=1)
    keep = lengths > 0
    lengths, is_nuc = lengths[keep], is_nuc[keep]

    seg_end = np.cumsum(lengths)
    seg_start = seg_end - lengths
    chain_end = np.cumsum(chain_length)
    chain_start = chain_end - chain_length

    return Template(seg_start, seg_end, is_nuc, chain_start, chain_end)
