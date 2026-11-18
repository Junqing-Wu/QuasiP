# Clamped–pinned von Kármán beam

Section 4.2, Fig. 14 and Table 3 of *Parallelized computation of quasi-periodic solutions for finite element problems: A Fourier series expansion-based shooting method* describe a clamped–pinned beam
with a linear spring and a periodic point load at midspan. The matching
SSMTool example is `examples/vonKarmanBeamIRs`. Its `build_model.m`
assembles the bare beam. Its `vonKarmanBeamIRs.m` adds the 37 N/mm spring
and sets the stiffness-proportional damping and 20 N midspan load.
The directory name contains `CC`; the boundary conditions implemented
in this example are clamped–pinned as specified in Section 4.2.

## Run in MATLAB

```matlab
nElements = 32;
[M,C,K,F,E,output] = build_vonKarmanBeamIRs(nElements);
```

The builder uses MATLAB sparse matrices and direct finite-element
assembly. Run `vonKarmanBeamIRs.m` to assemble the 32-element model
and print the first two undamped natural frequencies. No SSMTool
startup or add-on is needed.

`M*qdd + C*qd + K*q + F(q) = E*cos(theta)`, with
`theta = Omega*t`. The harmonic phase origin is arbitrary.

| Output | Definition |
|---|---|
| `M,C,K` | Reduced mass, damping, and stiffness. The spring is included in `K`, and `C=(2e-4/9)*K`. |
| `F` | Quadratic and cubic internal force handle. `F(Q)` takes one displacement state per column. `[f,J]=F(q)` also returns the nonlinear tangent for one state. |
| `E` | Midspan excitation vector including the 20 N amplitude. |
| `output` | Two rows that recover transverse displacement at one-quarter span and midspan. |

The axial field uses linear interpolation and the transverse field
uses cubic Hermite interpolation. Each node has axial displacement,
transverse displacement, and rotation. The left three DOFs and right
two translational DOFs are constrained. With 32 elements the model
has 94 free DOFs. Spring, force, and sensor locations are evaluated
at their physical positions using the same Hermite interpolation.

For a QuasiP callback that expects the right-hand side of the
second-order equation, use `E*cos(theta)-F(q)`.

## SSMTool-compatible nonlinear tensors

With MATLAB Tensor Toolbox on the path, request the same sparse tensor
layout as SSMTool's `fnl`:

```matlab
[M,C,K,F,E] = build_vonKarmanBeamIRs(32,'tensor');
f = double(ttv(F{1},{q,q},[2,3])) + double(ttv(F{2},{q,q,q},[2,3,4]));
```

`F{1}` is an `n x n x n` quadratic tensor and `F{2}` is an
`n x n x n x n` cubic tensor. Input indices are sorted and coefficients
of repeated monomials are summed. The QuasiP builder performs the
assembly locally. Tensor Toolbox provides the `sptensor` representation
and contraction. Omitting the format argument keeps the MATLAB-only
function-handle interface.

## Periodic response in QuasiP

Follow the [QuasiP installation instructions](../../README.md#installation), then make this example directory the MATLAB working directory.

```matlab
demo_BeamIRs_PO;
figure_BeamIRs_PO;
```

`demo_BeamIRs_PO.m` computes the full 94-DOF periodic branch with
MP/mVCF harmonic balance. The excitation is 20 N. COCO corrects the
linear complex response at 33 rad/s and continues in `k_OM1` toward
34.5 rad/s. The script uses harmonic order 7 and 32 phase samples.
Results are saved in `data/beamIRs_po_stab`. The figure script plots
quarter-span and midspan amplitudes and the largest Floquet multiplier
modulus. Solid and dashed sections show stable and unstable periodic
solutions. Fold, period-doubling (PD) and Neimark-Sacker (TR) labels
are displayed when detected.

Floquet analysis uses 256 steps per period. The demo sets the COCO event
tolerance to `1e-10` to resolve the small TR monitor values for this
model. With `store_all_psi=false`, the branch data store multipliers
at every point and the full transition history at TR points for restarts.
The `L` callback supplies sparse matrices. The MP stability routine
converts the matrix used for its rank check and inversion locally.

`BeamIRs_L` returns the constant matrices. `BeamIRs_F` applies the
harmonic load and subtracts the elastic nonlinear force. `BeamIRs_DFDX`,
`BeamIRs_DFDP`, and `BeamIRs_DFDTH` supply analytical derivatives.
The periodic solver uses the element force handle and its tangent.

## Quasi-periodic response after TR

After the periodic run, execute:

```matlab
demo_BeamIRs_QP_SHT;
demo_BeamIRs_QP_mVCF;
figure_BeamIRs_QP;
```

The SHT demo first calls `ode_TR2QuasiP` at the lower periodic TR and
corrects the seed with a fixed projection of the first sideband. It
saves the corrected seed in `data/beamIRs_qp_sht_seed`, releases the
projection constraint, and continues the torus branch to locate three
frequency events. The branch is saved in `data/beamIRs_qp_sht`.

The SHT discretization uses 1024 Newmark steps per excitation period,
32 transverse phases, and harmonic order 7 in the second angle. It
stores the full-state spectra for the mVCF starting guesses. The SHT
and mVCF demos run sequentially to limit RAM use. Quasi-periodic
stability analysis is disabled.

The mVCF demo corrects the SHT results at 33.859, 33.931, and
34.034 rad/s at fixed excitation frequencies. The mVCF harmonic orders
are (11,7), (13,7), and (11,7), with a 32-by-32 equation grid as in
the paper. Its sparse Jacobian preserves the assembled beam's local
coupling structure. The plots compare the quarter-span amplitude and
the second response frequency.



