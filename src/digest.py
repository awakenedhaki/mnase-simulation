from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple

import numpy as np

from src.chromatin import Template


@dataclass(frozen=True)
class DigestParams:
    """Kinetic parameters for MNase digestion.

    Rates are per unit time (arbitrary, e.g. minutes). Endonuclease rates are
    per bp; exonuclease rates are bp removed per fragment end.

    Attributes:
        k_linker (float): Endonuclease cut rate in linker DNA (per bp).
        k_nuc (float): Endonuclease cut rate inside nucleosome cores (per bp).
        k_exo_linker (float): Exonuclease trimming rate of a free end in linker DNA.
        k_exo_nuc (float): Exonuclease trimming rate of a free end inside a core.
        dt (float): Time step of the tau-leaping simulation.
        min_len (int): Fragments shorter than this are lost (acid-soluble).
    """

    k_linker: float = 1e-3
    k_nuc: float = 1e-5
    k_exo_linker: float = 1.0
    k_exo_nuc: float = 0.02
    dt: float = 0.1
    min_len: int = 20


def endonuclease_cuts(
    template: Template, cut_rate: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """Draws the cut sites produced during one time step.

    Cuts form a Poisson process whose intensity is piecewise constant over
    segments. A cut at position p breaks the bond between bp p-1 and bp p.

    Args:
        template (Template): The chromatin map.
        cut_rate (np.ndarray): Expected number of cuts per segment in this step.
        rng (np.random.Generator): Random number generator.

    Returns:
        np.ndarray: Sorted, unique cut positions.
    """
    cumulative = np.cumsum(cut_rate)
    n_cuts = rng.poisson(cumulative[-1])
    if n_cuts == 0:
        return np.empty(0, dtype=np.int64)
    seg = np.searchsorted(cumulative, rng.random(n_cuts) * cumulative[-1], side="right")
    cuts = rng.integers(template.seg_start[seg] + 1, template.seg_end[seg] + 1)
    return np.unique(cuts)


def apply_cuts(
    starts: np.ndarray, ends: np.ndarray, cuts: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Splits fragments at the given cut sites.

    Cuts that fall in DNA that has already been digested away, or exactly on an
    existing fragment end, are ignored. Because fragments are disjoint and each
    remaining cut lies strictly inside one of them, sorting the new starts and
    ends separately keeps them paired.

    Args:
        starts (np.ndarray): Sorted fragment starts.
        ends (np.ndarray): Sorted fragment ends (exclusive).
        cuts (np.ndarray): Cut positions.

    Returns:
        Tuple[np.ndarray, np.ndarray]: The new (starts, ends).
    """
    idx = np.searchsorted(starts, cuts, side="right") - 1
    clipped = np.clip(idx, 0, None)
    inside = (idx >= 0) & (cuts > starts[clipped]) & (cuts < ends[clipped])
    cuts = cuts[inside]
    return np.sort(np.concatenate([starts, cuts])), np.sort(np.concatenate([ends, cuts]))


def exonuclease_trim(
    template: Template,
    starts: np.ndarray,
    ends: np.ndarray,
    params: DigestParams,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray]:
    """Trims both ends of every fragment inward.

    Each end advances by a Poisson number of bp at the rate of the region it
    sits in, and never past the edge of that region in a single step, so ends
    stall at nucleosome borders before slowly nibbling into the core.

    Args:
        template (Template): The chromatin map.
        starts (np.ndarray): Sorted fragment starts.
        ends (np.ndarray): Sorted fragment ends (exclusive).
        params (DigestParams): Kinetic parameters.
        rng (np.random.Generator): Random number generator.

    Returns:
        Tuple[np.ndarray, np.ndarray]: The trimmed (starts, ends); fragments may
            now be empty or inverted and must be filtered by the caller.
    """
    rates = np.array([params.k_exo_linker, params.k_exo_nuc]) * params.dt

    is_nuc, _, seg_end = template.region_at(starts)
    step = rng.poisson(rates[is_nuc.astype(int)])
    starts = starts + np.minimum(step, seg_end - starts)

    is_nuc, seg_start, _ = template.region_at(ends - 1)
    step = rng.poisson(rates[is_nuc.astype(int)])
    ends = ends - np.minimum(step, ends - seg_start)

    return starts, ends


def simulate(
    template: Template,
    params: DigestParams,
    sample_times: Iterable[float],
    rng: Optional[np.random.Generator] = None,
) -> Dict[float, np.ndarray]:
    """Digests the template with MNase and samples fragment lengths over time.

    Every lane is an aliquot of the same digestion, so later samples are
    further-digested versions of earlier ones, as in a real time course.

    Args:
        template (Template): The chromatin map.
        params (DigestParams): Kinetic parameters.
        sample_times (Iterable[float]): Times at which to take an aliquot.
        rng (Optional[np.random.Generator]): Random number generator.

    Returns:
        Dict[float, np.ndarray]: Fragment lengths (bp) for each sample time.
    """
    rng = np.random.default_rng() if rng is None else rng
    sample_steps = {t: int(round(t / params.dt)) for t in sorted(sample_times)}
    cut_rate = (
        np.where(template.seg_is_nuc, params.k_nuc, params.k_linker)
        * template.seg_length
        * params.dt
    )

    starts, ends = template.chain_start.copy(), template.chain_end.copy()
    samples = {}
    for step in range(max(sample_steps.values()) + 1):
        for t, sample_step in sample_steps.items():
            if sample_step == step:
                samples[t] = ends - starts

        cuts = endonuclease_cuts(template, cut_rate, rng)
        starts, ends = apply_cuts(starts, ends, cuts)
        starts, ends = exonuclease_trim(template, starts, ends, params, rng)

        keep = ends - starts >= max(params.min_len, 1)
        starts, ends = starts[keep], ends[keep]

    return samples
