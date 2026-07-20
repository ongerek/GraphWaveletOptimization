# GraphWaveletOptimization — Statistically Optimized Filter Banks

**From autocorrelation-adapted compaction to precision-driven graph
wavelets.** Pure MATLAB / GNU Octave, no toolboxes required in Octave.
(also with an alternative working name `statfb`.)

[![CI](https://github.com/ongerek/GraphWaveletOptimization/actions/workflows/ci.yml/badge.svg)](https://github.com/ongerek/GraphWaveletOptimization/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Octave](https://img.shields.io/badge/GNU%20Octave-%E2%89%A5%206-blue)](https://octave.org)

`GraphWaveletOptimization` implements, in three layered modules, the
statistical design of critically sampled perfect-reconstruction (PR)
multiscale transforms:

1. **Module 1 — LP-optimal orthonormal filter banks.** Length-$N$
   paraunitary two-channel banks with $p$ vanishing moments maximizing
   energy compaction for a given input autocorrelation. Orthonormality
   is exact (half-band equalities), vanishing moments are structural
   ($F = B^pT$), the objective is linear, and the spectral factorization
   is a failure-proof cepstral method. Half-band error of delivered
   filters: ~1e-17.
2. **Module 2 — Wiener lifting.** Biorthogonal lifting with closed-form
   constrained-Wiener predict and Gouze (reconstruction-optimal) update;
   exact detail autocorrelation; Gauss–Markov performance bounds via
   precision-matrix conditioning. Includes the universality result that
   the moment-constrained two-tap design equals LeGall 5/3 for every
   WSS input.
3. **Module 3 — Graph lifting pyramids.** Statistics-adapted lifting on
   arbitrary weighted graphs: max-cut bipartition, precision
   conditional-mean / hop-restricted Wiener / $k$-nearest-neighbor
   Wiener / harmonic predictors, matrix Gouze update with masked closed
   form, exact covariance propagation, and Kron-reduction coarsening
   (Laplacian-preserving and equal to the exact marginal precision).
   Includes convex Laplacian-constrained graph learning. Reduces
   **exactly** to Module 2 on a path graph — verified to 1e-15.

## Repository layout

| Path | Contents |
|---|---|
| `src/` | all core functions (19 files) |
| `tests/` | verification harnesses: `run_tests` (Module 1), `run_tests2` (Module 2), `run_tests3` (Module 3) — 31 automated pass/fail checks in total; `run_all_tests` runs everything and errors on any failure |
| `demos/` | `demo_benchmark`, `demo_lifting_benchmark`, `demo_graph_benchmark` — the three benchmark studies |
| `paper/` | LaTeX manuscript (`main.tex`, `refs.bib`); `runme.m`, which regenerates every figure and table of the paper into `paper/paper_out/`; and `fig8_localization_sweep.m`, a standalone generator for the manuscript's support-size sweep figure (Fig. 8) |
| `.github/workflows/` | CI: runs the full test suite in GNU Octave on every push |

## Requirements

- **GNU Octave ≥ 6** — no extra packages (the LP uses Octave's built-in
  GLPK), *or*
- **MATLAB** — base MATLAB plus the Optimization Toolbox (`linprog`);
  the wrapper `src/lp_solve.m` selects the available solver
  automatically.

## Quick start

```matlab
% from the repository root
addpath src

% Module 1: adapt a 16-tap, 2-vanishing-moment orthonormal bank to AR(1)
r   = acf_from_psd(@(w) psd_ar([1 -0.95], 1-0.95^2, w));
out = design_compaction_fb(r, 16, 2);
G   = coding_gain_ortho(@(rr) out.h0, r, 4)     % 4-level coding gain, dB

% Module 2: Wiener lifting, 6 predict / 6 update taps
wl  = design_wiener_lifting(r, 6, 6);
Gb  = coding_gain_bior(wl.h0, wl.h1, wl.g0, wl.g1, r, 4)

% Module 3: 3-level pyramid on a graph GMRF (L: Laplacian, Sg: covariance)
pyr = graph_lifting_pyramid(L, Sg, 3, struct('predict','wiener','hops',1));
Gg  = tc_gain(pyr.T, pyr.S, Sg)

% Module 3, k-NN-restricted predictor (support-size sweep of Fig. 8):
% per odd vertex, the k nearest even vertices in weighted shortest-path
% distance
pyk = graph_lifting_pyramid(L, Sg, 3, struct('predict','wiener_knn','knn',6));
```

## Verification

Every mathematical claim in the paper marked "verified" has an
automated check:

```matlab
addpath src tests
run_all_tests    % 31 automated checks; errors if any fails
```

Checks include: spectral factorization and Daubechies / CDF 9/7 / LeGall
5/3 recovery against literature coefficients (~1e-15); perfect
reconstruction of every extracted and pyramid transform (~1e-15);
half-band certificates of adapted designs (~1e-17); the AR(1)
Markov-blanket collapse and the GMRF conditional-variance bound
(~1e-16); Kron reduction = marginal precision; and Monte-Carlo
agreement of the analytic detail variance. Several checks certify more
than one identity per line (e.g., off-diagonal sign and zero row sums
of the Kron-reduced Laplacian).

## Reproducing the paper

```matlab
cd paper
addpath ../src
runme                       % all figures (PNG+EPS) and LaTeX table
                            % fragments -> paper/paper_out/  (a few minutes)
fig8_localization_sweep     % Fig. 8 only (support-size sweep, sensor net)
```

### A note on cross-platform reproducibility

All classical (Modules 1–2), path-graph, and graph-learning numbers
reproduce digit-for-digit across MATLAB and Octave builds. The
sensor-network gains (Table VII and Fig. 8) may vary by roughly
±0.01 dB across platforms: the greedy max-cut bipartition is seeded by
the sign pattern of a Laplacian eigenvector, and the sign returned by
`eig` is LAPACK-build-dependent, which can select a different (equally
valid) local optimum of the cut. The graph itself, its GMRF, and the
KLT bound are identical everywhere; only bipartition-downstream third
decimals move.

## Citing

If you use this code, please cite the accompanying manuscript,
*"Statistically Optimized Lifting Filter Banks: From
Autocorrelation-Adapted Compaction to Precision-Driven Graph Wavelets"*
(see `CITATION.cff`; a preprint link will be added upon posting).

## License

MIT — see [LICENSE](LICENSE).
