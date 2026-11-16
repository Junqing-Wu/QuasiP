# QuasiP

QuasiP is a MATLAB toolbox built on [COCO](https://sourceforge.net/projects/cocotools/) for computing and continuing periodic and quasi-periodic solutions of nonlinear ordinary differential equations. It supports first- and second-order systems, including autonomous systems and systems under periodic or quasi-periodic excitation.

A response with $d$ base frequencies is represented by a $d$-dimensional invariant torus. The excitation prescribes $e$ of those frequencies ($0\leq e\leq d$), and QuasiP solves for the remaining $d-e$ frequencies together with the torus. The same framework covers periodic solutions with $d=1$ and tori with more than two base frequencies.

## Capabilities

- **Continuation and bifurcation analysis.** Follow periodic or quasi-periodic solution branches as model parameters or excitation frequencies vary. COCO provides continuation, monitoring, event localization, and restart. Periodic-orbit analysis includes fold, period-doubling, and Neimark–Sacker events. A saved Neimark–Sacker point can initialize a two-frequency torus.
- **Torus doubling.** Independent period-multiplication factors $P_i$ extend the torus domain in each angular coordinate. This allows doubled tori and their continuations to be computed within the same framework.
- **Unknown frequencies and phase conditions.** QuasiP applies an orthogonal phase condition for each unknown base frequency using the continuation iterates. The user need not supply a reference torus.
- **Response and stability.** Reconstruct amplitudes, trajectories, and frequency spectra from the computed solutions. Optional post-processing uses Floquet multipliers for periodic orbits and Lyapunov exponents for quasi-periodic tori.

## Numerical methods

| Method | Formulation | Typical use |
|---|---|---|
| **m-VCF** (multi-step variable-coefficient formulation) | Uses a multi-period parametrization and discretizes the PDEs for the full invariant torus. Supports harmonic balance (HB), collocation (CO), and finite differences (FD), with alternating evaluation in the coefficient and sample domains. | Direct torus computation and detailed studies of discretization, amplitudes, and spectra. |
| **FSE-SHT** (Fourier series expansion-based shooting) | Uses a multi-trajectory parametrization. Integrates trajectories from a torus cross-section and imposes phase-shifted matching at the end of a cycle. Supports hyper-time (HT) and multi-frequency (MF) shooting formulations. | Large models, including finite-element systems. The cross-section formulation reduces the number of shooting unknowns, and trajectory integrations can run in parallel. |

Both methods generate nonlinear algebraic equations for COCO. The toolbox follows a **Starter–Generator–Processor** workflow. The Starter builds the initial problem, the Generator supplies the discretized equations and phase conditions, and the Processor reconstructs responses and evaluates stability.

## Installation

1. Download [COCO](https://sourceforge.net/projects/cocotools/) separately. In MATLAB, go to the downloaded COCO directory and run `startup`.
2. In the same MATLAB session, go to the QuasiP root directory and run `install`.

Run each example from its own directory. COCO normally stores results under `data/<run_name>/` in the working directory. Parallel trajectory integration and some example-specific routines require the corresponding MATLAB toolboxes. See each example's README for its run order and dependencies.

## Starter functions

| Entry point | Purpose |
|---|---|
| `ode_isol2QuasiP` | Construct a continuation problem from model callbacks, an initial Fourier representation, frequencies, and parameters. |
| `ode_QuasiP2QuasiP` | Restart from a saved QuasiP solution. |
| `ode_TR2QuasiP` | Predict a two-frequency torus from a saved periodic-orbit Neimark–Sacker point, labeled `TR`. |
| `ode_BP2QuasiP` | Start a secondary branch using the restart tangent at a saved COCO branch point, labeled `BP`. |

## Examples

| Example | What it demonstrates |
|---|---|
| [Duffing oscillator](example/0001_Duffing/README.md) | Quasi-periodic forcing, comparison of discretizations, Lyapunov stability, torus doubling, and continuation with multiple base frequencies. |
| [Self-excited airfoil with a store](example/0002_Airfoil/README.md) | Autonomous periodic and two-frequency responses, initialization at a torus bifurcation, and comparison of m-VCF and FSE-SHT. |
| [Rotor–stator rubbing](example/0003_RSrubbing/README.md) | Periodic and quasi-periodic response curves, stability, and orbit comparisons at prescribed rotor speeds. |
| [Clamped–pinned von Kármán beam](example/0004_CC_vonKarman_beam/README.md) | Full-order finite-element assembly, periodic Floquet analysis, and an FSE-SHT/m-VCF comparison for the two-frequency response after a Neimark–Sacker bifurcation. |

## Documentation and references

The toolbox manuscript, *QuasiP: A COCO-Based Toolbox for Continuation of Quasi-periodic Solutions in Nonlinear ODEs* by Junqing Wu, Jun Jiang, and Mingwu Li, describes the problem formulation, numerical methods, implementation, and examples. See Sections 2–4 for the formulation and toolbox structure and Section 5 for the numerical examples.

- Junqing Wu, Ling Hong, Mingwu Li, and Jun Jiang. *General multi-steps variable-coefficient formulation for computing quasi-periodic solutions with multiple base frequencies*. **Acta Mechanica Sinica**, 42, 425381 (2026).
- Junqing Wu, Ling Hong, Mingwu Li, and Jun Jiang. *Parallelized computation of quasi-periodic solutions for finite element problems: A Fourier series expansion-based shooting method*. **Mechanical Systems and Signal Processing**, 259, 114779 (2026). 
- Junqing Wu, Ling Hong, Yeyin Xu, and Jun Jiang. *An efficient and robust approach for continuation and bifurcation analysis of quasi-periodic solutions by multi-harmonic balance method*. **Journal of Sound and Vibration**, 601, 118943 (2025). 
- Junqing Wu, Ling Hong, and Jun Jiang. *A comparative study on multi- and variable coefficient harmonic balance methods for quasi-periodic solutions*. **Mechanical Systems and Signal Processing**, 187, 109929 (2023).
- Junqing Wu, Ling Hong, and Jun Jiang. *A robust and efficient stability analysis of periodic solutions based on harmonic balance method and Floquet-Hill formulation*. **Mechanical Systems and Signal Processing**, 173, 109057 (2022).
- Mingwu Li. *Tor: a toolbox for the continuation of two-dimensional tori in autonomous systems and non-autonomous systems with periodic forcing*. arXiv:2012.13256 (2020). [Preprint](https://arxiv.org/abs/2012.13256) · [Tor source code](https://github.com/mingwu-li/torus_collocation)
- Robert Fiedler, Hartmut Hetzler, and Simon Bäuerle. *Efficient numerical calculation of Lyapunov-exponents and stability assessment for quasi-periodic motions in nonlinear systems*. **Nonlinear Dynamics**, 112, 8299–8327 (2024).
- Harry Dankowicz and Frank Schilder. *Recipes for Continuation*. SIAM (2013).

## Acknowledgments

We would like to thank Hongming Liang and Feng Wang (Southern University of Science and Technology) for valuable discussions during the development of QuasiP. 
