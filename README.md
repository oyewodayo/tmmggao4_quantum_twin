# Quantum twin of TmMgGaO₄ — CERN QTI Hackathon

Working implementation of the full challenge pipeline: warm-up AFM state
preparation (Phase 1), a Rydberg-array "quantum twin" of the frustrated
triangular magnet TmMgGaO₄ reproducing Fig. 1 of Leclerc et al.
(arXiv:2603.20372) (Phase 2a), and an open-exploration gentle quench /
thermalisation study (Phase 2b).

All three notebooks in `notebooks/` are **pre-executed** with real
outputs from this codebase (small system sizes, chosen to run on a
laptop CPU in under two minutes each) — open them to see the actual
figures, or re-run everything yourself.

## Layout

```
src/material_mapping.py   -- all physics: device/register builders,
                              the TmMgGaO4 <-> Rydberg Hamiltonian
                              mapping (Eq. 1/3/4 of the paper), pulse
                              sequence builders for all three phases,
                              and observable helpers.
notebooks/
  01_phase1_afm_prep.ipynb              -- Scholl et al. AFM warm-up
  02_phase2a_material_twin.ipynb        -- TmMgGaO4 magnetisation curve
  03_phase2b_quench_thermalization.ipynb -- gentle quench (open Phase 2b)
requirements.txt
```

## Quickstart

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
jupyter lab notebooks/
```

Everything runs on CPU at the sizes shipped here (N=9 for Phase 1,
N=25 for Phase 2a/2b). This is intentional so you can `git clone` and
get a green run immediately, and so the physics/API usage is validated
before you burn GPU time.

## Scaling up on a GPU

This is where the real hackathon points live. The knobs to turn, all
exposed as plain function arguments in `src/material_mapping.py`:

| Knob | Where | Small demo | Paper-scale target |
|---|---|---|---|
| Register size | `triangular_rhombus_register(l_bulk, r1, buffer_rows)` | `l_bulk=3, buffer_rows=1` (N=25) | `l_bulk=3, buffer_rows=2` (N=49, paper's smallest); `l_bulk=6/9` for N=100/169 |
| Bond dimension | `emu_mps.MPSConfig(max_bond_dim=...)` | 12, 24 | 128, 256, ... — **always run at least 2 values and check convergence**, this is an explicit evaluation criterion |
| Sweep duration | `t_rise, t_sweep, t_fall` in the sequence builders | ~1.5 μs | longer sweeps improve adiabaticity for bigger systems, at the cost of more MPS timesteps |
| Post-quench hold time | `t_hold_over_J1` in notebook 3 | 3 | 10+ to see the classical-frontier bond-dimension blowup (Fig. 4/Ext. Dat. Fig. 6 of the paper) |

`emu_mps.MPSConfig(num_gpus_to_use=1)` (or more) puts the TDVP solver on
GPU — see the `emu-mps` docs for multi-GPU options.

## What's deliberately left as an exercise

* **Experimental data overlay** (Phase 2a, Fig. 1e comparison): we do
  not fabricate digitised values of the paper's AC-susceptibility curve.
  Notebook 2 has a ready-made loader (`material_data.csv`) — digitise
  the published figure yourself (e.g. WebPlotDigitizer) and drop the
  file in to auto-overlay it.
* **QMC-thermal comparison** (Phase 2b, Fig. 4c-style dashed lines):
  needs a stochastic-series-expansion sampler for the Rydberg
  Hamiltonian (Sandvik, Phys. Rev. E 68, 056701 (2003)) — a good
  "further work" callout for your writeup, out of scope for a
  hackathon-time implementation.
* **AFM structure factor** for Phase 1 (currently we only compute
  staggered magnetisation): one line via
  `pulser.backend.CorrelationMatrix`, following the same pattern used
  for `C1^zz` in notebook 3.

## Key references

- L. Leclerc et al., *One-to-one quantum simulation of the
  low-dimensional frustrated quantum magnet TmMgGaO4 with 256 qubits*,
  arXiv:2603.20372 (2026). — Phase 2 anchor, reproduce its Fig. 1.
- P. Scholl et al., *Programmable quantum simulation of 2D
  antiferromagnets with hundreds of Rydberg atoms*, Nature 595, 233
  (2021), arXiv:2012.12268. — Phase 1 anchor.
- J. Vovrosh et al., arXiv:2511.19340 and arXiv:2511.20388 — the
  classical-frontier motivation for Phase 2b.
- H. Silvério et al., *Pulser*, Quantum 6, 629 (2022).
