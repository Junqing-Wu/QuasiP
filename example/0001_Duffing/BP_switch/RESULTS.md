# Executed BP switching test

Tested with MATLAB R2022b and the installed COCO/QuasiP on 2026-11-14.
The original Duffing demos, their data, and toolbox functions were unchanged.

## Detected branch points

The Tor11 seed from `qp514_tor11` was corrected and continued in the
P=(2,2) equations. The run `tor11_P22` covered omega1 from 1.8 to 3
with 129 branch-data points and detected:

| BP label | omega1 | Saved BP vector norm |
|---|---:|---:|
| 37 | 1.92689674939579 | 1 |
| 109 | 2.83087842035248 | 1 |

These values are numerical locations for the stated discretization,
not a claim of that many digits of physical accuracy. Their proximity
to the previously observed Tor22 endpoints supports testing the
Tor11-to-Tor22 transition as BP switching in the enlarged representation.

## Branch identity check

The diagnostic is the Euclidean norm of saved displacement Fourier
coefficients with noninteger physical wave numbers. It is evaluated
for both coordinate directions as well as their union. The doubled-box
Tor11 control contains only integer wave numbers.

| Run | Points | omega1 coverage | Maximum noninteger coefficient norm |
|---|---:|---|---:|
| `tor11_P22` | 129 | [1.8, 3] | 0 |
| `switch_BP37` | 63 | [1.861628, 2.108802] | 0 |
| `switch_BP109` | 63 | [2.625281, 2.989414] | 0 |
| `direct_BP37` | 28 | [1.926897, 1.929962] | 0.231648 |
| `direct_BP109` | 28 | [2.829736, 2.830879] | 0.233209 |

For each direct run, the maximum directional norms in directions 1
and 2 equal the listed union norm (to the displayed precision). Both
directions therefore acquire doubled-period content. The direct runs
are local Tor22 branches, whereas the original-mode controls remain
Tor11. The diagnostic uses the stored spectra; zero means no such modes
were retained in those spectra, rather than an exact-arithmetic claim.

The lower direct branch develops toward increasing omega1 and the
upper direct branch toward decreasing omega1. Both continuation signs
were attempted. The tiny frequency excursions near the initial BP,
below 4e-7, should not be interpreted as a resolved additional fold.

## Why the original call returned Tor11

`QuasiP_read_solution(...,'BP')` already retrieves `lsol.v` and passes
it to the constructed problem as `sol.t0`. This vector represents the
transverse near-null direction of the bordered system.
The original `ode_BP2QuasiP` also requests COCO `branch='switch'`.
In this COCO version, `CurveSegment.init_TS` then computes a transverse
direction relative to the supplied tangent, which can select the
primary Tor11 branch again. The two original-mode controls did so.

The example-local wrapper `ode_BP2QuasiP_direct.m` calls the original
function and then sets:

```matlab
prob = coco_set(prob,'cont','branch','no switch');
```

This tells COCO to refine and use the supplied secondary vector. Both
direct-mode runs acquired substantial half-integer content without
an integration-based Tor22 seed. The wrapper is deliberately local;
the original toolbox implementation has not been replaced.

## Scope and reproduction

Run the scripts in README.md, including the optional control, to reproduce the comparison.
`figure_5_1_4_BP` draws the amplitude comparison and the noninteger
coefficient-norm comparison.mat`.
COCO logs and saved solutions are in the five subdirectories of `data`.

This establishes a local numerical branch switch at both endpoints in
the tested P=(2,2) discretization. It does not establish full Tor22
coverage in one run, resolution convergence of BP locations, or
Lyapunov stability of the switched solutions. Stability processing
was disabled during this algebraic test. Only isolated boundary points
are BPs; the entire unstable Tor11 interval is not a branch point.
