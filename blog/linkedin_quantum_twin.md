Exactly two weeks ago, our team sat down for a CERN Quantum Technology Initiative hackathon with a deceptively simple question: can a quantum computer simulate a real crystal, atom for atom — not a toy model, an actual material you could put under a microscope?

The challenge: build a "quantum twin" of TmMgGaO₄, a frustrated magnet where spins sit on a triangular lattice and simply can't all satisfy their neighbors at once (imagine seating three mutually feuding guests around one small table — someone's always unhappy). That frustration is exactly what makes these materials scientifically interesting, and exactly what makes them hard to simulate classically.

Neutral-atom Rydberg hardware implements the material's governing equations natively — no digital approximation needed — which is what makes it a genuine "quantum twin" rather than a generic model. The challenge was emulation-only (no QPU access required), so we built and validated everything on Pasqal's own emulators — exact simulation for small systems, tensor-network emulation (emu-mps) on real Pasqal Cloud jobs for larger ones. Concretely, we:

→ Warmed up by reproducing a landmark 2D antiferromagnet result (Scholl et al., Nature 2021), then pushed it further with a real 25-atom run on Pasqal's cloud backend
→ Mapped the real TmMgGaO₄ Hamiltonian onto a 49-atom Rydberg array and reproduced its published magnetization curve, digitizing the real experimental data ourselves and overlaying it directly against our simulation
→ Pushed into quench dynamics and found the exact point where classical simulation software (tensor networks) starts to break down under growing entanglement — the frontier where a quantum device stops being a nice-to-have and becomes the only tool that works

The biggest lesson here wasn't a physics result — it was what actually reproducing a published paper takes. A figure in a paper is the destination, not the road: it hides a hundred small decisions the text never spells out — undocumented pulse timings, unstated conventions, device parameters that don't quite match your setup. Rebuilding it from scratch, and being honest about exactly where and why our version diverges, taught us more about the physics than reading the paper ever did.

Full technical writeup and code linked below.

https://medium.com/@oyewodayo/quantum-twin-of-a-frustrated-magnet-simulating-a-real-frustrated-magnet-tmmggao%E2%82%84-thulium-4c14e2c33dc8

#QuantumComputing #CERN #Physics #QuantumSimulation #NeutralAtoms #CondensedMatter
