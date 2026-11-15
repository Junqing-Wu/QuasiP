# Duffing: Tor11 embedded in the Tor22 problem and BP switching

This independent experiment tests whether the algebraic branch points
of an embedded Tor11 can be used by the existing `ode_BP2QuasiP`
to initialize a doubled torus. It does not modify the original Duffing
scripts or their saved continuation runs.

## Run order

Follow the [QuasiP installation instructions](../../../README.md#installation), then make this example directory the MATLAB working directory.

```matlab
demo_5_1_4_BP_a;
demo_5_1_4_BP_b;
figure_5_1_4_BP;
```

To reproduce the original-mode comparison, run `demo_5_1_4_BP_control`
after `demo_5_1_4_BP_a` and before plotting. Each demo writes its run data
to the `data` directory under `BP_switch`.

The prerequisite is `../data/qp514_tor11`, produced by
`demo_5_1_4a`. That saved branch supplies only the starting Tor11 Fourier
coefficients. Its first saved point is at omega1=1.8. No integrated Tor22
guess, perturbation copied from an existing Tor22, or segmented Tor22
continuation is used to initialize branch switching.

## What each script does

1. `demo_5_1_4_BP_a.m` reconstructs the Tor11 initial guess using its physical
   wave numbers, embeds it in the same P=(2,2), MP/mVCF, HB/HB setup used
   by `demo_5_1_4b`, and continues it with COCO BP detection enabled.
   The run is `tor11_P22`. The frequency ratio omega2=omega1/sqrt(2) and
   model parameters `(zeta,gamma,f1,f2)=(0.1,0.2,5,5)` remain fixed.
2. `demo_5_1_4_BP_control.m` calls the existing `ode_BP2QuasiP` for each saved
   BP, reattaches the frequency-ratio equation, and attempts both branch
   directions. Runs are named `switch_BP<LAB>`. This is the optional control and can be run after the scan.
3. `demo_5_1_4_BP_b.m` calls the toolbox ode_BP2QuasiP to follow the saved BP null vector.
   The short runs are named `direct_BP<LAB>` (12 steps per direction).
4. `figure_5_1_4_BP.m` compares amplitudes and the norm of Fourier
   coefficients with noninteger physical wave numbers. It prints the numerical diagnostics and draws two figures.

The input HB orders are (7,7), S=(32,32), S_wave=(256,256), and
P=(2,2). QuasiP expands the period-dependent basis internally; inspect
the printed equation table for the actual number of coefficients and
grid points. Do not double the physical `k_wave` values when embedding
the original torus: changing P already provides the enlarged domain.

## Interpreting the test

Tor11 is a subspace of the P=(2,2) equations. Integer physical wave
numbers alone still describe Tor11. Nonzero half-integer content in
both directions shows that the switched solution has left that subspace
and has doubled periods in both directions. Inspect the coefficient
norm along the switched branch, rather than only its first BP point,
where the doubled-torus amplitude is zero.

An entire unstable interval is not a set of branch points. Candidate BPs
are isolated intersections of solution branches where the discretized
Jacobian becomes singular. A change in a saved Lyapunov stability flag
does not by itself establish an algebraic BP or guarantee a switchable
Tor22 branch.

Lyapunov post-processing is disabled for this initial algebraic test to
avoid evaluating long QR iterations at every scan point. Consequently
these new runs do not carry physical stable/unstable labels. The scan
does not establish convergence of stability exponents or of the
infinite-dimensional bifurcation location. Inspect the saved BP
frequencies, the noninteger coefficient norms, and the continuation logs
before interpreting a successful switch.

## COCO tangent handling in this installation

`QuasiP_read_solution(...,'BP')` retrieves `lsol.v`, the null vector of
the bordered system, and `ode_construct_QuasiP` supplies it as `t0`.
The previous implementation additionally selected `branch='switch'`.
In the installed COCO implementation (`CurveSegment.init_TS`), this mode
computes another transverse direction relative to the supplied tangent.
The repaired toolbox function selects `branch='no switch'`, so
COCO refines and follows the supplied secondary vector directly.
The control demo explicitly restores the previous COCO mode.

Numerical findings from the original 0004 comparison are recorded in
`RESULTS.md`.
