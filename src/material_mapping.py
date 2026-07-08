"""
material_mapping.py
====================

Core toolkit for the CERN QTI hackathon challenge "Quantum twin of a
frustrated magnet: TmMgGaO4 on a neutral-atom processor".

This module implements the physics that ties the two phases of the
challenge together:

  * Phase 1  -- a plain square-lattice Rydberg array realising the
                textbook 2D antiferromagnetic Ising model (Scholl et al.,
                Nature 595, 233 (2021), arXiv:2012.12268).

  * Phase 2  -- a triangular Rydberg array whose native Hamiltonian is
                mapped onto the microscopic model of TmMgGaO4, following
                Leclerc et al., "One-to-one quantum simulation of the
                low-dimensional frustrated quantum magnet TmMgGaO4 with
                256 qubits", arXiv:2603.20372 (2026).

Everything here is *emulation-only*: sequences are built with Pulser and
meant to be run on pulser_simulation.QutipBackend (small systems, exact)
or emu_mps.MPSBackend (larger systems, approximate -- always check
bond-dimension convergence, see notebooks/02 and 03).

Conventions / units
--------------------
All energies (J1, Delta_x, Delta_z, Omega, delta) are expressed as
*angular* frequencies in rad/us, matching Pulser's internal convention.
"h-bar * J1" in the paper is what we call `J1` here (i.e. we always work
with hbar = 1, exactly as the paper does after writing H/hbar).

Hamiltonian recap (Leclerc et al., Eqs. 1, 3, 4)
-------------------------------------------------
Material (triangular-lattice transverse-field Ising model):

    H_TMGO / hbar = J1 * sum_<ij> sz_i sz_j + J2 * sum_<<ij>> sz_i sz_j
                     + sum_i (Delta_x * sx_i - Delta_z * sz_i)

QPU (Rydberg Hamiltonian, ground state |g> = up, Rydberg state |r> = down):

    H_QPU / hbar = sum_{i<j} U_ij n_i n_j + (Omega(t)/2) sum_i sx_i
                   - delta(t) sum_i n_i,          n_i = (1 - sz_i) / 2

    U_ij = C6 / r_ij^6

Identifying nearest-neighbour interaction with J1:

    J1 = U_1 / 4 = C6 / (4 * r1^6)                      (r1 = NN spacing)

and the control-parameter mapping (Eq. 4, up to O(next-nearest-neighbour)
corrections H_diff that we neglect here, exactly as most of the main text
figures do in their "quasi-classical" / leading-order treatment):

    Delta_x(t) = Omega(t) / 2
    Delta_z(t) = (delta_U - delta(t)) / 2,   delta_U = (1/2) sum_j U_1j / N

For a bulk (translation-invariant) site, delta_U ~= (1/2) * sum over all
neighbour shells of U_ij, which for a large lattice is dominated by the
6 nearest neighbours: delta_U ~= 3 * U_1 = 12 * J1. We use this bulk
approximation (also stated explicitly in the paper's Methods: "This
interaction profile is constant in the bulk of the triangular lattice").
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from pulser import Register, Sequence, Pulse
from pulser.waveforms import RampWaveform, ConstantWaveform, CompositeWaveform
from pulser.devices import VirtualDevice
from pulser.channels import Rydberg


# ---------------------------------------------------------------------------
# 1. Device construction
# ---------------------------------------------------------------------------

def build_device(rydberg_level: int = 75, min_atom_distance: float = 4.0) -> VirtualDevice:
    """Build a VirtualDevice matching one of the paper's Orion Beta QPUs.

    Pasqal's Orion Beta devices used in the paper:
        FM1 (France):     |75 S_1/2, m_J=1/2>   -- used for most data
        FM2/FC1:           |60 S_1/2, m_J=1/2>

    Passing `rydberg_level=75` reproduces C6/(2*pi*hbar) ~= 1949 GHz.um^6,
    matching the value quoted in the paper's Methods for FM1.

    We use a VirtualDevice (no hardware calibration limits on Rabi
    frequency, detuning range or atom number) since this challenge is
    emulation-only and we want to match the *physical* Hamiltonian of the
    paper exactly, not a specific hardware guardrail.
    """
    return VirtualDevice(
        name=f"orion_beta_like_n{rydberg_level}",
        dimensions=2,
        rydberg_level=rydberg_level,
        min_atom_distance=min_atom_distance,
        channel_objects=(Rydberg.Global(max_abs_detuning=None, max_amp=None),),
    )


def j1_from_device(device: VirtualDevice, r1: float) -> float:
    """hbar*J1 (rad/us) = C6 / (4 * r1^6), Eq. 4 of Leclerc et al."""
    return device.interaction_coeff / (4.0 * r1**6)


# ---------------------------------------------------------------------------
# 2. Registers
# ---------------------------------------------------------------------------

def square_register(rows: int, cols: int, spacing: float, prefix: str = "q") -> Register:
    """Plain square lattice, for Phase 1 (Scholl et al. AFM warm-up)."""
    return Register.rectangle(rows, cols, spacing=spacing, prefix=prefix)


def triangular_rhombus_register(
    l_bulk: int, r1: float, buffer_rows: int = 2, prefix: str = "q"
) -> tuple[Register, np.ndarray]:
    """Triangular-lattice rhombus register, following the paper's Ext. Dat.
    Fig. 3 / Methods prescription.

    The paper arranges atoms on an N = L x L rhombus, with L a multiple of
    3 plus `buffer_rows` extra rows on each side of the bulk so that
    boundary effects don't contaminate bulk observables (L = 7, 10, 13, 16
    for bulk sizes 3, 6, 9, 12).

    Parameters
    ----------
    l_bulk : int
        Linear size of the bulk region of interest. Must be a multiple of
        3 for commensurability with the 1/3-ordered phase.
    r1 : float
        Nearest-neighbour spacing (um).
    buffer_rows : int
        Extra rows of atoms added on *each* side to reduce edge effects
        (paper default: 2).

    Returns
    -------
    register : pulser.Register
    is_bulk : np.ndarray[bool]
        Boolean mask (same order as register.qubits) marking which atoms
        belong to the (L_bulk x L_bulk) bulk region used for observables,
        as opposed to the boundary buffer.
    """
    if l_bulk % 3 != 0:
        raise ValueError(
            f"l_bulk={l_bulk} must be a multiple of 3 for commensurability "
            "with the 1/3-filling order (see Leclerc et al., Methods)."
        )
    L = l_bulk + 2 * buffer_rows

    coords = []
    is_bulk = []
    lo, hi = buffer_rows, buffer_rows + l_bulk  # bulk index range [lo, hi)
    for row in range(L):
        for col in range(L):
            x = col * r1 + row * (r1 / 2.0)
            y = row * r1 * np.sqrt(3) / 2.0
            coords.append((x, y))
            is_bulk.append(lo <= row < hi and lo <= col < hi)
    coords = np.array(coords)
    coords -= coords.mean(axis=0)
    reg = Register.from_coordinates(coords, prefix=prefix)
    return reg, np.array(is_bulk, dtype=bool)


def nearest_neighbour_bonds(register: Register, r1: float, tol: float = 0.15):
    """Return list of (i, j) qubit-index pairs at distance ~r1 (bulk NN bonds)."""
    names = list(register.qubits.keys())
    coords = np.array([register.qubits[n] for n in names], dtype=float)
    bonds = []
    n = len(names)
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(coords[i] - coords[j])
            if abs(d - r1) < tol * r1:
                bonds.append((i, j))
    return bonds


# ---------------------------------------------------------------------------
# 3. Material <-> QPU parameter mapping (Eq. 4)
# ---------------------------------------------------------------------------

@dataclass
class MaterialParams:
    """TmMgGaO4 microscopic model parameters, from Leclerc et al. (citing
    Li et al. PRX 10, 011007 (2020) and H. Li et al. Nat. Commun. 11, 1111
    (2020) for the original characterisation)."""

    J2_over_J1: float = 0.05          # next-nearest-neighbour coupling
    Delta_x_over_J1: float = 1.08     # intrinsic transverse (crystal) field
    delta_z_per_field: float = 1.543  # Delta_z/J1 per Tesla of mu0*H


MATERIAL = MaterialParams()


def delta_U_bulk(J1: float, n_nn: int = 6) -> float:
    """Bulk site-averaged interaction sum delta_U = (1/2) sum_j U_1j,
    dominated by the `n_nn` nearest neighbours (6 on a triangular
    lattice): delta_U ~= (n_nn/2) * U_1 = 2 * n_nn * J1 (since U_1 = 4*J1).
    """
    return (n_nn / 2.0) * (4.0 * J1)


def qpu_controls_from_material(
    Delta_x_over_J1: float, Delta_z_over_J1: float, J1: float
) -> tuple[float, float]:
    """Material (Delta_x/J1, Delta_z/J1) -> QPU controls (Omega, delta), rad/us.

    Inverting Eq. 4:
        Delta_x = Omega / 2                => Omega = 2 * Delta_x
        Delta_z = (delta_U - delta) / 2    => delta = delta_U - 2*Delta_z
    """
    Omega = 2.0 * Delta_x_over_J1 * J1
    dU = delta_U_bulk(J1)
    delta = dU - 2.0 * Delta_z_over_J1 * J1
    return Omega, delta


def material_from_qpu_controls(Omega: float, delta: float, J1: float) -> tuple[float, float]:
    """Inverse of the above: QPU controls -> (Delta_x/J1, Delta_z/J1)."""
    Delta_x_over_J1 = (Omega / 2.0) / J1
    dU = delta_U_bulk(J1)
    Delta_z_over_J1 = (dU - delta) / (2.0 * J1)
    return Delta_x_over_J1, Delta_z_over_J1


def mu0H_from_Delta_z(Delta_z_over_J1: float, params: MaterialParams = MATERIAL) -> float:
    """Convert Delta_z/J1 to the physical applied field mu0*H (Tesla),
    using the paper's quoted Delta_z/J1 ~= 1.543 * mu0*H(T)."""
    return Delta_z_over_J1 / params.delta_z_per_field


# ---------------------------------------------------------------------------
# 4. Analytic classical reference (zero transverse field, Delta_x -> 0)
# ---------------------------------------------------------------------------

def classical_Mz(Delta_z_over_J1: np.ndarray, J2_over_J1: float = MATERIAL.J2_over_J1):
    """Classical (Delta_x = 0) magnetisation, from the energy comparison in
    the main text: E_up...up/N = 3*J1 - Delta_z + O(J2), Mz=1, versus
    E_1/3/N = -J1 - Delta_z/3 + O(J2), Mz=1/3. The crossover sits at
    Delta_z/J1 ~= 6 (quoted in the text for J2 ~ 0.05*J1). This is a sharp
    step in the strict classical limit -- useful as a sanity-check
    reference curve, *not* a fit to data.
    """
    Delta_z_over_J1 = np.asarray(Delta_z_over_J1, dtype=float)
    Dz_c = 6.0  # crossover quoted in the paper for J2/J1 ~ 0.05
    return np.where(Delta_z_over_J1 < Dz_c, 1.0 / 3.0, 1.0)


# ---------------------------------------------------------------------------
# 5. Pulse sequence builders
# ---------------------------------------------------------------------------

def afm_prep_sequence(
    register: Register,
    device: VirtualDevice,
    Omega_max: float,
    delta_start: float,
    delta_end: float,
    t_rise: float,
    t_sweep: float,
    t_fall: float,
) -> Sequence:
    """Phase-1 style quasi-adiabatic AFM-preparation sequence (Scholl et
    al. 2021 protocol): ramp Omega up while sweeping delta from large
    negative to positive/zero, then ramp Omega back down to freeze the
    z-basis populations before readout.

    Timing mirrors Ext. Dat. Fig. 4a of Leclerc et al.: Omega ramps up
    over `t_rise`, both pulses run together for `t_sweep`, Omega ramps
    down over `t_fall` while delta keeps sweeping to `delta_end`.
    """
    seq = Sequence(register, device)
    seq.declare_channel("ryd", "rydberg_global")

    amp = CompositeWaveform(
        RampWaveform(t_rise, 0.0, Omega_max),
        ConstantWaveform(t_sweep, Omega_max),
        RampWaveform(t_fall, Omega_max, 0.0),
    )
    det = RampWaveform(t_rise + t_sweep + t_fall, delta_start, delta_end)
    seq.add(Pulse(amp, det, 0.0), "ryd")
    return seq


def quasi_adiabatic_material_sweep(
    register: Register,
    device: VirtualDevice,
    J1: float,
    Omega_max_over_J1: float,
    Delta_z_start_over_J1: float,
    Delta_z_end_over_J1: float,
    t_rise: float,
    t_sweep: float,
    t_fall: float,
) -> Sequence:
    """Phase-2a sequence: single continuous quasi-adiabatic sweep across
    the material's antiferromagnetic quantum phase transition, built
    directly in *material* units and converted to QPU controls via the
    Eq. 4 mapping. Sampling <Mz(t)> along this one sweep (and converting
    t -> Delta_z(t)/J1 via the mapping) reproduces a full magnetisation
    curve like Fig. 1e / Ext. Dat. Fig. 4a, in a single simulation.
    """
    Omega_max = Omega_max_over_J1 * J1
    _, delta_start = qpu_controls_from_material(0.0, Delta_z_start_over_J1, J1)
    _, delta_end = qpu_controls_from_material(0.0, Delta_z_end_over_J1, J1)

    seq = Sequence(register, device)
    seq.declare_channel("ryd", "rydberg_global")
    amp = CompositeWaveform(
        RampWaveform(t_rise, 0.0, Omega_max),
        ConstantWaveform(t_sweep, Omega_max),
        RampWaveform(t_fall, Omega_max, 0.0),
    )
    det = RampWaveform(t_rise + t_sweep + t_fall, delta_start, delta_end)
    seq.add(Pulse(amp, det, 0.0), "ryd")
    return seq


def gentle_quench_sequence(
    register: Register,
    device: VirtualDevice,
    J1: float,
    Omega_max_over_J1: float,
    Delta_z_prep_over_J1: float,
    Delta_z_quench_over_J1: float,
    t_prep_rise: float,
    t_prep_sweep: float,
    t_hold: float,
) -> Sequence:
    """Phase-2b sequence: adiabatically prepare the ground state at
    Delta_z_prep/J1 (paramagnetic side, following the paper's post-quench
    protocol which starts from |up...up>), then abruptly (square-pulse)
    switch the longitudinal field to Delta_z_quench/J1 and hold, so the
    resulting unitary evolution can be probed for thermalisation of local
    observables (cf. Fig. 4 of Leclerc et al.).
    """
    Omega_max = Omega_max_over_J1 * J1
    _, delta_prep = qpu_controls_from_material(0.0, Delta_z_prep_over_J1, J1)
    _, delta_quench = qpu_controls_from_material(0.0, Delta_z_quench_over_J1, J1)

    seq = Sequence(register, device)
    seq.declare_channel("ryd", "rydberg_global")

    # Preparation: bring the system close to |g...g> = |up...up> ground state
    # of a strongly paramagnetic point, then quench.
    amp_prep = CompositeWaveform(
        RampWaveform(t_prep_rise, 0.0, Omega_max),
        ConstantWaveform(t_prep_sweep, Omega_max),
    )
    det_prep = RampWaveform(t_prep_rise + t_prep_sweep, delta_prep, delta_prep)
    seq.add(Pulse(amp_prep, det_prep, 0.0), "ryd")

    # Square-pulse quench: instantaneous change of the longitudinal field,
    # amplitude held constant (post-quench Hamiltonian sampled at fixed
    # Omega, delta -- as in Ext. Dat. Fig. 4b).
    seq.add(
        Pulse(ConstantWaveform(t_hold, Omega_max), ConstantWaveform(t_hold, delta_quench), 0.0),
        "ryd",
    )
    return seq


# ---------------------------------------------------------------------------
# 6. Observable helpers
# ---------------------------------------------------------------------------

def sz_from_occupation(n: np.ndarray) -> np.ndarray:
    """Convert Rydberg-state occupation n_i = (1 - sz_i)/2 to sz_i."""
    return 1.0 - 2.0 * np.asarray(n)


def bulk_magnetisation(n: np.ndarray, is_bulk: Optional[np.ndarray] = None) -> float:
    """Average sz over the bulk region (or all sites if is_bulk is None)."""
    sz = sz_from_occupation(n)
    if is_bulk is not None:
        sz = sz[is_bulk]
    return float(np.mean(sz))
