# QuasiP: a COCO-based toolbox for continuation of quasi-periodic solutions with multiple base frequencies in quasi-periodically forced ODEs

We present a COCO-based toolbox for the continuation of quasi-periodic solutions with $d$ base frequencies in nonlinear ODEs subjected to quasi-periodic excitation with $e$ base frequencies, where $d \geq e \geq 0$. To address the computation of $d$-dimensional invariant tori corresponding to quasi-periodic solutions in the phase space, we define two types of hyper-time domains for torus parametrization, namely multi-periodic parameterization and multi-trajectory parameterization, thereby deriving the boundary value problems (BVPs) of multi-variable functions governed by transformed PDEs and ODEs, respectively. For the first parameterization, we propose a single-step constant-coefficient formulation (\textit{s}-CCF) and a multi-step variable-coefficient formulation (\textit{m}-VCF) from the perspectives of full discretization and semi-discretization, respectively. These formulations extend and unify the nonlinear algebraic equations (NAEs) required by the harmonic balance method, collocation method, and finite difference method for solving invariant tori of arbitrary dimensions. For the second parameterization, we also propose a single-step variable-coefficient formulation (\textit{s}-VCF) to construct the NAEs for the collocation method, mapping method, and shooting method in solving invariant tori of arbitrary dimensions. Additionally, we propose a simple and robust phase condition to fix the phase shift of the torus induced by unknown base frequencies.

*Currently, only \textit{m}-VCF, PC and stability analysis are available.

Please refer to [1-4] for more details about \textit{s}-CCF, \textit{m}-VCF, \textit{s}-VCF and phase condition. For stability analysis, refer to [5-6].

This toolbox is developed in the platform COCO (continuation core). For more details about COCO, please refer to [7-8]. 

When you use this toolbox and want to cite it, you may use the following



## References

[1] Wu JQ, Hong L, Jiang J. A comparative study on multi- and variable coefficient harmonic balance methods for quasi-periodic solutions [J]. Mechanical Systems and Signal Processing, 2023,
187:109929.

[2] Wu JQ, Hong L, Li MW, Jiang J. General multi-steps variable-coefficient formulation for computing quasi-periodic solutions with multiple base frequencies [J]. Acta Mechanica Sinica, 2026, 42(1):
425381.

[3]. Wu JQ, Hong L, Li MW, Jiang J. Parallelized computation of quasi-periodic solutions for finite element problems: A Fourier serie sexpansion-based shooting method [J]. Mechanical Systems and Signal Processing, MSSP25-3979R1, in Second Review.

[4] Wu JQ, Hong L, Xu YY, Jiang J. An efficient and robust approach for continuation and bifurcation
analysis of quasi-periodic solutions by multi-harmonic balance method [J]. Journal of Sound and
Vibration, 2025, 601:118943.

[5] Wu JQ, Hong L, Jiang J. A robust and efficient stability analysis of periodic solutions based on harmonic balance method and Floquet-Hill formulation [J]. Mechanical Systems and Signal Processing,
2022, 173:109057.

[6] FIEDLER R, HETZLER H, B?UERLE S. Efficient numerical calculation of Lyapunov-exponents and stability assessment for quasi-periodic motions in nonlinear systems[J]. Nonlinear Dynamics, 2024, 112: 8299-8327.

[7] https://sourceforge.net/projects/cocotools/

[8] Dankowicz, H., & Schilder, F. (2013). Recipes for continuation. Society for Industrial and Applied Mathematics.
