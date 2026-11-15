# Duffing oscillator example

## Overview

This example follows Section 5.1, **"A SDOF Duffing oscillator system,"**
of the QuasiP article. The governing equation is

```text
xddot + 2*zeta*xdot + x = -gamma*x^3 + sum_i f_i*sin(omega_i*t).
```

For the main two-frequency examples, `zeta=0.1`, `gamma=0.2`, `f1=f2=5`,
and `omega2=omega1/sqrt(2)`. Each script specifies its own initial frequency
and parameter interval. The experiments demonstrate equivalent ODE forms,
initialization and restart, numerical-method comparison, Lyapunov stability,
torus doubling, and extension to more base frequencies.

Figure numbers below refer to the supplied article. Local script numbering
is different from the article figure numbering.

## 1. Setup

Follow the [QuasiP installation instructions](../../README.md#installation), then make this example directory the MATLAB working directory.

Keep this example folder as the working directory for computation and
plotting. With default COCO settings, results are saved in
`data/<run_name>/`. Results generated in another working directory are
not automatically available here. The MATLAB functionality used by
`findpeaks` and the parallel stability routines must also be available.

## 2. Recommended execution order

| Group | Scripts in order | Article content |
|---|---|---|
| A | `demo_5_1_1a` -> `demo_5_1_1b` -> `figure_5_1_1` | Eqs. (37)-(39), Fig. 6: equivalent formulations and restart |
| B | `demo_5_1_2` -> `figure_5_1_2` | Fig. 7: comparison of five numerical configurations |
| C | `demo_5_1_3` -> `figure_5_1_3` | Fig. 8: response branch and Lyapunov stability |
| D/E | `demo_5_1_4a` -> `demo_5_1_5` -> `figure_5_1_5` -> `demo_5_1_4b` -> `demo_5_1_4c` -> `figure_5_1_4` | Figs. 9-10: torus doubling, phase portraits, and spectra |
| F | `demo_5_1_6a`, `demo_5_1_6b`, `demo_5_1_6c` | Fig. 11: two-, three-, and four-base-frequency solutions |

Groups A-C are independent. Group F is optional and each script includes
its own plotting. For group D/E, the order is intentionally not strictly
numerical: `demo_5_1_5` prepares initial solutions needed by `4b/4c`.

```matlab
% Equivalent formulations and restart (Fig. 6)
demo_5_1_1a;
demo_5_1_1b;
figure_5_1_1;

% Numerical-method comparison (Fig. 7)
demo_5_1_2;
figure_5_1_2;

% Lyapunov stability (Fig. 8)
demo_5_1_3;
figure_5_1_3;

% Torus branches and representative solutions (Figs. 9-10)
demo_5_1_5;
figure_5_1_5;    % Plot immediately while integration arrays remain available

demo_5_1_4a;
demo_5_1_4b;
demo_5_1_4c;
figure_5_1_4;
```

Several demos call `clearvars` or `close all`. Export figures before
starting the next group if they must be preserved.

## 3. Experiments and saved runs

### A. Equivalent formulations and restart - Fig. 6

`demo_5_1_1a` implements the three equivalent forms in Eqs. (37)-(39):
the scalar second-order equation, the explicit first-order state equation,
and a matrix first-order equation. Initial guesses come from numerical
integration. `ode_isol2QuasiP` constructs continuation in `omega1`.

`demo_5_1_1b` reads saved `Om1=3` events and uses `ode_QuasiP2QuasiP`
to restart in `gamma`. Run `1a` first and ensure those events were saved.
Outputs are `qp511a_sys1/2/3` and `qp511b_sys1/2/3`.

### B. Numerical-method comparison - Fig. 7

`demo_5_1_2` compares five configurations: m-VCF with HB/FD/CO choices,
HT shooting with Newmark, and MF shooting with `ode45`.
`figure_5_1_2` compares displacement and velocity time histories.
These are zero-dimensional COCO solves at a selected parameter point,
not five full continuation branches.
Outputs are `qp512_meth1` through `qp512_meth5`.

### C. Stability analysis - Fig. 8

`demo_5_1_3` generates `qp513_stability`.
`figure_5_1_3` displays response amplitudes and Lyapunov information.
The choices `Omt=1` and `Omt=2` evaluate stability using different
base-frequency directions. The branch columns `Stab` and `Ly` store
the stability labels and exponents, respectively.

### D. Torus-doubling branches - Fig. 9

Tor11, Tor22, and Tor44 use period factors `P=(1,1)`, `(2,2)`, and `(4,4)`.
They all remain two-base-frequency computations.

| Script | Required saved input | Output runs |
|---|---|---|
| `demo_5_1_4a` | None; its own initial guess | `qp514_tor11` |
| `demo_5_1_4b` | `qp515_tor22` for its final section | `qp514_tor22_a`, `qp514_tor22_b`, `qp514_tor22_bb` |
| `demo_5_1_4c` | `qp515_tor44` | `qp514_tor44_a`, `qp514_tor44_b` |
| `figure_5_1_4` | All the `qp514_*` runs above | Combined branch plots |

Run `demo_5_1_5` before executing the full `4b` and before `4c`.
The table describes the segmented `4b` currently in this folder.
If a single-run version replaces it, update the Tor22 inputs of
`figure_5_1_4` to read that complete branch once instead of the three
segmented runs.

### E. Phase portraits and spectra - Fig. 10

`demo_5_1_5` computes representative Tor11, Tor22, and Tor44 solutions
at initial `omega1=2.85`, `2.75`, and `2.71`, respectively.
The current script produces six single-point runs:

```text
qp515_tor11       qp515_tor11_b
qp515_tor22       qp515_tor22_b
qp515_tor44       qp515_tor44_b
```

The `_b` runs provide the additional time-domain sampling used by the
plotting script. They are generated by the current demo and are not
missing prerequisites.

Run `figure_5_1_5` immediately after `demo_5_1_5`, in the same MATLAB
workspace. Besides saved runs, the figure uses direct-integration arrays
such as `x_tor11`, `X_tor11`, and `ff_tor11`, with their Tor22/Tor44
counterparts. It compares numerical integration with reconstructed phase
portraits, Poincare sections, and spectra, as discussed in Section 5.1.

### F. More base frequencies - Fig. 11

`demo_5_1_6a/b/c` illustrate `d=e=2/3/4` at `omega1=2.85`.
Each script constructs its own initial guess, performs a single-point
solve, and plots phase portraits, sections, and spectra.
Outputs are `qp516_e2`, `qp516_e3`, and `qp516_e4`.
The four-frequency dense reconstruction can require substantial memory.

## 4. Current code notes

- `demo_5_1_6a` still lists `k_OM3` in its COCO parameter list although
  `d=e=2`. Remove this nonexistent parameter name before running it.
- `demo_5_1_6c` lists only three frequency names. Include `k_OM4` if
  all four known frequencies should be explicitly reported.
- `figure_5_1_1` selects `labOm13(4:5)`. If a new run has fewer than
  five `Om1=3` events, adapt this selection to the available labels.
- Inspect each `coco_log.txt` and actual branch coverage. Requesting an
  interval does not guarantee that the entire interval was reached.
- As stated in Section 5.1, Lyapunov exponents assess stability but do
  not by themselves identify torus-doubling bifurcations. These examples
  use selected period factors and spectral initial guesses to compute
  the different torus families.

Checked against the current scripts and the supplied QuasiP article,
Section 5.1, on 14 November 2026. This is an execution guide, not a claim
that every current demo and figure has been rerun successfully.
