# Airfoil/store example — QuasiP Section 5.2

This example studies the autonomous three-coordinate airfoil/store model
in Section 5.2, Eqs. (40)–(44), of the QuasiP article. The coordinates are
relative plunge q1, wing pitch q2, and store pitch q3. The equations have
the form `Mq*qddot + Cq*qdot + Kq(Q)*q + [0;0;10*q3^3] = 0`.
Q is the continuation parameter. There are no prescribed excitation
frequencies (`e=0`); response frequencies are solved as unknowns.

## Setup and recommended order

Follow the [QuasiP installation instructions](../../README.md#installation), then make this example directory the MATLAB working directory.

```matlab
demo_5_2_p;
demo_5_2_qp_mVCF;
demo_5_2_qp_SHT;

close all;
fig_5_2_1;
close all;
fig_5_2_2;
close all;
fig_5_2_3;
```

The figure scripts reuse figure numbers. Inspect or save each set before
closing its windows and drawing the next set. COCO normally writes runs
to `data/<run_name>/` under the working directory. Existing complete runs
can be plotted without repeating continuation.

## Computation scripts and dependencies

| Script | Required run | Output and purpose |
|---|---|---|
| `demo_5_2_p.m` | None | Integrates at Q=6, estimates a periodic Fourier guess, and continues `qp_d1` in Q over the requested interval [5,10]. Uses MP/mVCF with HB order 9. |
| `demo_5_2_qp_mVCF.m` | `qp_d1` with TR labels | Starts from the smallest saved TR label with epsilon=0.1. Computes `qp_d2_mvf_a` near the lower boundary and the longer run `qp_d2_mvf_b`, using MP/mVCF and HB orders (9,7). |
| `demo_5_2_qp_SHT.m` | `qp_d1` with TR labels | Starts from the largest saved TR label with epsilon=0.1. Computes `qp_d2_sht_a`, using MT/SHT in the MF domain, `ode15s`, and orders (9,5). |
| `demo_5_2_3.m` | `qp_d1` with TR labels | Constructs `prob2` from a TR point for inspection. It has no active `coco` call and does not generate a branch or a figure. |

The article reports Hopf onset near Q=4.7440 and torus bifurcations near
Q=7.0762 and Q=9.5035. The periodic demo starts from integration at Q=6;
it does not continue the equilibrium branch from the Hopf point.
The current scripts select TR points using the minimum/maximum saved
label. Check their actual Q values before identifying them as the lower
and upper physical boundaries. Label order alone is not a physical test.

## Plotting scripts

| Script | Required runs | Contents |
|---|---|---|
| `fig_5_2_1.m` | `qp_d1`, `qp_d2_mvf_a`, `qp_d2_mvf_b`, `qp_d2_sht_a` | Three displacement-amplitude diagrams against Q, comparing periodic and quasi-periodic branches and both discretizations; corresponds to Fig. 12. |
| `fig_5_2_2.m` | Same four runs | Response frequencies against Q and the frequency–frequency diagram; corresponds to Fig. 13. |
| `fig_5_2_3.m` | `qp_d1` with both TR points and saved `Psi`, `HT_x`, `HT_u` | At each TR point overlays the periodic orbit and predicted tori with epsilon=0.05, 0.10, 0.20 in `(q1,qdot1,q2)`; corresponds to Fig. 14. |

`fig_5_2_3` calls `Airfoil_TR_predict.m`, an example-local copy of the
prediction routines. Plotting is implemented here, outside
`ode_TR2QuasiP`. The colored tori are uncorrected initial guesses from the
critical Floquet mode, not converged quasi-periodic solutions. This
figure needs only the periodic run, so it may be drawn immediately after
`demo_5_2_p`.

## Model callbacks and numerical interpretation

`Airfoil1_F.m` supplies the nonlinear force, `Airfoil1_DFDX.m` its state
derivative, and `Airfoil1_L.m` the linear model matrices. The periodic
problem has one unknown frequency and one phase condition; the
quasi-periodic problem has two of each.

Periodic stability is evaluated by Floquet analysis. Quasi-periodic
stability uses saved Lyapunov estimates; the active QP demos use
`S=32`, `L1=256`, `N_Ly=500` (mVCF: both Omt directions; SHT: Omt=1).
Near a boundary, check convergence of near-zero exponents before making
a stability claim. A requested continuation interval or point budget
does not guarantee coverage: inspect `data/<run_name>/coco_log.txt`.

A sharp turn near a TR point should be checked against residuals, step
sizes, frequencies, and nonzero second-frequency Fourier content. A
small QP amplitude may converge back to the embedded periodic solution;
a visually smooth curve or a successful corrector alone cannot exclude
that possibility.
