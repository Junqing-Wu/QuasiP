# Rotor-stator rubbing example: event solutions

This example implements the two-degree-of-freedom rotor–stator rubbing
system in Section 5.2 of Wu et al. (2025),
[*An efficient and robust approach for continuation and bifurcation analysis of quasi-periodic solutions by multi-harmonic balance method*](https://doi.org/10.1016/j.jsv.2025.118943),
Journal of Sound and Vibration, 601, 118943.

Follow the [QuasiP installation instructions](../../README.md#installation), then make this example directory the MATLAB working directory.

```matlab
demo_TR2QP;
figure_RS_FRC;
figure_RS_Lyapunov;
figure_RS_events;
```

`demo_TR2QP` registers the following events on both the periodic branch
`qp1` and the quasi-periodic branch `qp2`:

| Event name | Prescribed excitation frequency |
|---|---|
| `Om1=0.30` | `k_OM1=0.30` |
| `Om1=0.35` | `k_OM1=0.35` |
| `Om1=0.40` | `k_OM1=0.40` |
| `Om1=0.45` | `k_OM1=0.45` |

The glue constraint keeps the excitation frequency equal to the model
rotation parameter `OM`. Events are localized by COCO, rather than
selected as the nearest output point. Existing runs without these event
labels must be recomputed to obtain them. Actual availability depends on
whether continuation reaches each value; inspect `coco_log.txt`.

The plotting script follows the existing example style: separate `figure`
windows, direct `plot` calls, LaTeX axis labels, and `plot_setting_paper()`.
It plots only rotor-plane orbits `(x,y)`, grouped by excitation frequency.

The plot explicitly identifies unstable periodic orbits, stable
quasi-periodic tori, and unstable quasi-periodic tori. Unstable periodic
orbits are black dashed curves; stable QP tori are red point clouds;
unstable QP tori are blue point clouds. Stable periodic orbits, if present,
are black solid curves. Inconsistent/unknown stability is labeled
undetermined. Classification uses the saved Stab column, without filtering
out unstable solutions. All labels at the same event are included, so folded branches may
produce more than one solution at the same frequency. A missing-class
warning is issued when no unstable PO, stable QP, or unstable QP is saved
at a requested event. Registering an event does not generate an untracked
branch or guarantee all three classes exist at every requested frequency.

Figure titles identify frequency; legends identify solution type, stability,
and LAB. Plot data
are reconstructed from each saved coefficient vector; the plotting script
does not change saved solutions, rerun continuation, or automatically
export figures. No numerical-integration trace is presented as validation
of an unstable continued solution.


## Frequency-response curves

Run `figure_RS_FRC` after computing both branches. Figure 1 plots
`max(abs(x))` versus `k_OM1`; Figure 2 plots `max(abs(y))` versus `k_OM1`.
Black curves represent periodic orbits and red curves represent
quasi-periodic tori. Solid segments are stable, dashed segments are
unstable, and dotted segments have undetermined/inconsistent stability.
Labels come from the saved `Stab` column; exponents are not recomputed.
Continuation point order is preserved so folds remain correctly connected.
## Evolution of stability exponents

Run `figure_RS_Lyapunov` to plot the saved exponents without recomputing
continuation or stability. Figure 3 shows the real parts of the four
periodic-orbit Floquet exponents; Figure 4 shows the four QP Lyapunov
exponents for each saved Omt direction; Figure 5 compares the largest
state exponents. Zero is the stability reference. The additional 0.001
line is the current QP classification threshold, not an exact zero.
The horizontal axis defaults to k_OM1; set use_continuation_number=true
in the script for evolution along branch-data point order. The spectra
are ranked pointwise and do not track a single eigenmode through crossings.
## Shared XY plots

figure_RS_events overlays all periodic and quasi-periodic event solutions
at the same k_OM1 in one XY figure (101-104). Legends identify stability
and LAB. No spectrum figures are created.


