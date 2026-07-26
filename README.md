# Quantum twin of TmMgGaO₄ — CERN QTI Hackathon

Working implementation of the full challenge pipeline: warm-up AFM state
preparation (Phase 1), a Rydberg-array "quantum twin" of the frustrated
triangular magnet TmMgGaO₄ reproducing Fig. 1 of Leclerc et al.
(arXiv:2603.20372) (Phase 2a), and an open-exploration gentle quench /
thermalisation study that finds the classical-emulator breakdown point
(Phase 2b).

All three notebooks in `notebooks/` are **pre-executed** with real
outputs from this codebase — open them to see the actual figures. Some
sections were run locally on CPU in under two minutes; others were run
for real on **Pasqal Cloud** (see below) and their outputs are cached
results from an actual quantum-emulator job, not something you get for
free by just cloning the repo.

## Layout

Each notebook is standalone: it contains the device/register builders,
the TmMgGaO4 <-> Rydberg Hamiltonian mapping, pulse sequence builders,
and observable helpers it needs.

```
notebooks/
  01_phase1_afm_prep.ipynb              -- Scholl et al. AFM warm-up
  02_phase2a_material_twin.ipynb        -- TmMgGaO4 magnetisation curve
  03_phase2b_quench_thermalization.ipynb -- gentle quench + classical-frontier sweep
  material_data.csv                     -- digitised Fig. 1e AC-susceptibility overlay data
blog/
  medium_quantum_twin.md, linkedin_quantum_twin.md  -- writeup drafts
requirements.txt
credentials.example.yaml                -- template; copy to credentials.yaml and fill in
```

## What's local vs. what needs Pasqal Cloud

Not all of this runs for free on a laptop. Each notebook mixes a local,
credential-free section with a Pasqal Cloud section:

| Notebook | Local (CPU, no account needed) | Pasqal Cloud (needs real credentials) |
|---|---|---|
| `01_phase1_afm_prep` | Sections 1-4: exact `QutipBackendV2` emulation, `N=9` (3x3), plus an optional guarded `N=16` exact check | Section 5: `N=25` (5x5) on EmuMPS, `chi=50/100` |
| `02_phase2a_material_twin` | none — the whole notebook is cloud-only | Whole notebook: full `N=49` (7x7) register on EmuMPS, `chi=50` main scan / `chi=100` convergence check, analysed over the central 5x5=25-site bulk |
| `03_phase2b_quench_thermalization` | Part A: exact-diagonalisation thermal reference + local `emu_mps` quench, `N=9` | Part B: the paper's `N=49` register, quench sweep over `chi in {50,100,200,300}` — this is the section that finds the classical bond-dimension divergence |

The local sections run in well under two minutes and are a good
first sanity check. The cloud sections are where the headline results
(the Fig. 1-style magnetisation curve, and the classical-frontier
divergence) come from, and they were genuinely submitted to and run on
Pasqal's EmuMPS backend — their outputs in the notebooks are real job
results, cached so you don't need an account just to *read* the paper.

## Quickstart

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
jupyter lab notebooks/
```

To re-run the local sections, that's all you need. To re-run (or
extend) the Pasqal Cloud sections:

```bash
cp credentials.example.yaml credentials.yaml   # fill in your real account details
```

`credentials.yaml` is gitignored — never commit real credentials.
Each notebook has a `SUBMIT_NEW_JOBS` (or `SUBMIT_SWEEP`) safety switch,
off by default; flip it once to submit fresh jobs, wait for them to
show `DONE` on the Pasqal dashboard, then flip it back and run the
retrieval cells. Batch IDs are cached in the `*.json` files at the repo
root so a notebook never resubmits jobs it already has results for.

## Scaling further

The knobs that matter for pushing past what's shipped here:

| Knob | Where | Notes |
|---|---|---|
| Register size | `n_side` (01), fixed `L=7` (02), `l_bulk`/`buffer_rows` via `triangular_rhombus_register` (03 Part A) | 02 and 03-Part-B are already at the paper's `N=49`; growing further means a bigger `L`/`l_bulk` and a matching credentials-backed cloud run |
| Bond dimension | `max_bond_dim` / `chi` args to `EmulationConfig` | 02 checks `chi=50` vs `100`; 03 Part B sweeps `{50,100,200,300}` — **always compare at least two values and check convergence**, this is an explicit evaluation criterion |
| Sweep / hold duration | `t_rise, t_sweep, t_fall` (01, 02), `t_hold_over_J1` (03) | Longer sweeps improve adiabaticity; longer holds are what exposes the bond-dimension blowup in 03 Part B |

`EmulationConfig`/cloud submissions don't expose a GPU flag directly —
the compute happens on Pasqal's infrastructure once you submit.

## What's deliberately left as an exercise

* **QMC-thermal comparison** (Phase 2b Part B, Fig. 4c-style dashed
  lines): needs a stochastic-series-expansion sampler for the Rydberg
  Hamiltonian (Sandvik, Phys. Rev. E 68, 056701 (2003)) — a good
  "further work" callout, out of scope for a hackathon-time
  implementation.
* **AFM structure factor** for Phase 1 (currently we only compute
  staggered magnetisation): one line via
  `pulser.backend.CorrelationMatrix`, following the same pattern used
  for `C1^zz` in notebook 3.
* **Field-dependence of the classical breakdown** (Phase 2b Part B):
  the divergence-vs-field sweep does not cleanly peak near the
  expected transition (`Delta_z/J1 ~ 3.9`); notebook 3 documents this
  honestly rather than papering over it, with finite-size recurrence
  as the leading suspect.

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
