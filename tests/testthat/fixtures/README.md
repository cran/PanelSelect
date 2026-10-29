# Numerical fixtures and provenance

## 1. Independent likelihood and score fixtures

`likelihoods.R` contains small deterministic reference panels for the seven panel
models, with both signs of nonzero correlations where relevant. It stores the
parameters, data, integration orders, log likelihood, and natural score.

The original values were generated with independent base-R calculations.
The included `../reference/math.R` contains the corresponding likelihood and
score references. They construct a normal GH rule from a symmetric Jacobi
eigenproblem, evaluate the bivariate CDF using an adaptive Plackett integral,
derive observation scores separately, and integrate probabilities and scores
directly. Small moderate-parameter cases were checked against adaptive singleton
integrals and mathematical identities.

Fixed-rule cases use 7 outer nodes and 13 inner nodes. Refined cases use 31 outer
nodes and 25 inner nodes. Fixed comparisons use matching rules and orientation;
adaptive comparisons allow explicitly declared quadrature approximation error.
Changing assumptions or reference values requires mathematical review. The
implementation under test is not a source of expected values.

## 2. High-precision reference fixtures

The three Rmpfr fixtures are readable R text loaded with `dget()`. Decimal
literals retain 17 significant digits and round-trip exactly to the stored
doubles. Their metadata records the mathematical method, checks, precision,
software versions, generation date, field meanings and independence statement.

Values were computed from mathematical definitions with Rmpfr 1.1.2 and
MPFR 4.2.2 under R 4.5.3 on 2026-10-11; PanelSelect was never loaded. Inputs
are the exact stored IEEE doubles, converted with `mpfr(x, bits)`. Precision
comparisons and independent checks preceded writing, and reloading verified
identical doubles. Rmpfr is not required to run the package tests: tests read
the frozen values and do not regenerate them.

| Fixture | Contents | Precision and independent checks |
|---|---|---|
| `normal-tail.R` | log Φ and derivatives of order 0–5, the Mills ratio and x + M at 41 points; 36 deep-tail log-ratios | 1024 bits; 512 bits, closed forms and central differences |
| `bvn-tail.R` | log P(X < a, Y < b; r) and its first partials at 53 points | 384 bits; 256 bits, swapped ordering and Plackett's identity |
| `linear-marginal.R` | Linear outcome marginal, gradient and conditional person-effect moments for 3 cases, 12 panels | 512 bits; 256 bits, dense Cholesky and central differences |

### `normal-tail.R`

- `points`: `x` contains -1e4, -1e3, -300, -200, -100, -60, -50, -40, -35, -31,
  the switch points -30, -25 and -15 each with ±1e-6 offsets, -4 with ±1e-9 and
  ±1e-3 offsets, and further values from -29 to 8. The fields are
  `value` = F0 = log Φ(x), `first`, ..., `fifth` = F1, ..., F5 with
  Fk = d^k log Φ(x)/dx^k, `mills` = M = F1, and `remainder` = x + M
  evaluated in high precision (not x + fl(M)).
- `log_ratio`: 36 pairs with `x0` in {-30.5, -35, -50, -100, -1e3, -1e4},
  `x1` = fl(x0 + d) for d in {-0.5, -1e-3, -1e-9, 1e-9, 1e-3, 0.4}, and
  `value` = log Φ(x1) - log Φ(x0). Every argument is below -30.
- Method: log Φ from MPFR's correctly rounded `erfc`; M = exp(log φ - log Φ);
  the Leibniz recurrence M^(n+1) = -x M^(n) - n M^(n-1) - Σ C(n,k) M^(k) M^(n-k);
  Fk = M^(k-1). Checks: 512 against 1024 bits (largest relative difference
  2e-108); independently derived closed forms in (x, M), for example
  F4 = (3x - x³)M + (4 - 7x²)M² - 12xM³ - 6M⁴ (3e-270); and second-order central
  differences of log Φ alone with step 2^-60 at 1536 bits (4e-35, the truncation
  level).
- Implementation mapping: `normal_log_derivatives(x, 3L)` returns
  `value`, `first`, `second` and `third`; `normal_mills(x)` is `mills`;
  `pln_inner_mills(x)` returns `ratio` = `mills` and `curvature` = -`second`;
  `inner_taylor_mills_derivatives(z)` returns c(`first`, ..., `fifth`), the Mills
  ratio and its derivatives of order 0–4. The tail log-ratio inside
  `log_bvn_tail` is `log_ratio$value`.
- Coverage: deep-tail derivatives (orders ≤ 3 and orders 4–5),
  continuity across -4, the `normal_mills` NA guard (values at -40 and 0), and
  the tail test in `test-local-models.R`.

### `bvn-tail.R`

- `points`: 45 `tail` points with P < 1e-9 by design, using r in {±0.9999,
  ±0.99, ±0.9, ±0.5, ±0.1}, thresholds down to -40, mixed signs, both argument
  orders and positive-r thresholds near b = ra, where the conditional factor is
  far from one. Eight threshold points solve log(P/1e-9) = -1e-6
  (`threshold_below`) or +1e-6 (`threshold_above`) for
  (b, r) = (-2, 0.5), (-2, 0.9999), (1.5, -0.5) and (0, -0.9). This margin is
  far larger than pbivnorm's relative error near 1e-9.
- Fields: `a`, `b`, `r`, `log_p` = log P, the first partials
  `dlog_a` = φ(a)Φ((b - ra)/s)/P, `dlog_b` = φ(b)Φ((a - rb)/s)/P and
  `dlog_r` = φ2(a, b; r)/P, where s = sqrt((1 - r)(1 + r)), `kappa` =
  |a dlog_a| + |b dlog_b| + |r dlog_r|, and `region`.
- Method: adaptive Gauss–Legendre sums (nodes solved by Newton's method in
  mpfr) of the conditional integral ∫ φ(z)Φ((b - rz)/s) dz over z < a, with Φ
  from MPFR `erfc`. Breakpoints are placed geometrically around the peak of the
  log-concave integrand, and tangent-line bounds keep each discarded tail at
  most 1.2e-53 of P. Checks: 384 against 256 bits with different rule pairs and
  tolerances (largest relative difference in P 2.3e-45), the swapped ordering
  (5.5e-44), and Plackett's identity integrated from r = -1,
  P = max(0, Φ(a) - Φ(-b)) + (2π)^-1 ∫ exp(-(a+b)²/(8 sin²(t/2)) -
  (a-b)²/(8 cos²(t/2))) dt over 0 < t < acos(-r), which has no cancellation
  (1.9e-44).
- Conditioning: a double evaluation that perturbs its inputs by a relative ε
  changes log P by about ε·kappa. Here 41 of the 53 points have ε·kappa ≤ 1e-12
  (ε = 2.2e-16). Seven points with r ≤ -0.9 exceed 1e-11, up to 6e-6
  (kappa = 2.7e10) at (-38, 5, -0.9999), where log P ≈ -2.7e6. An absolute
  log-error test should therefore select points by `kappa` or allow a
  conditioning term.
- Implementation mapping: `log_bvn_value(a, b, r)` is `log_p`. The partials in
  `probit_re_probit_re_rows` (shared with Probit–Probit) and
  `prediction_bvn` (`d_A`,
  `d_signed_B`, `d_r`) are `dlog_a`, `dlog_b` and `dlog_r`.
- Coverage: batched tail integration (absolute log error on
  well-conditioned points), continuity across the 1e-9 switch, and argument
  symmetry.

### `linear-marginal.R`

- `cases`: `moderate` (σ = 1, λ = 0.5; m = 1, 2, 5, a panel with 2 of 4 rows
  selected, and an unselected panel), `intermediate` (σ = 1e-2, λ = 1;
  m = 50, 3), and `cancellation` (σ = 1e-3, λ = 2; m = 1, 2, 10, 50, 200 with
  residual means between 1.4 and 6.3).
- Each case stores `data` (`id`, `y` with NA where `d` = 0, `d`, `x`, `w`),
  `theta` (`outcome.(Intercept)`, `outcome.x`, `sd_outcome` = λ,
  `sd_error` = σ), `count` (selected rows per panel), and per panel in id order
  `value`, `gradient`, `mean`, `variance`, `scale`, `mean_gradient`,
  `variance_gradient` and `scale_gradient`. Gradients are n × 4 matrices whose
  columns are named by `theta`. Panels with m = 0 have value 0, zero gradients,
  mean 0, variance 1 and scale 1.
- Data lie on dyadic grids (x on 2^-20, y and the intercept on 2^-30, the slope
  on 2^-10), so the residuals y - Xβ are exact in double for any summation
  order. The fixture therefore tests the marginal formula, not residual
  rounding.
- Method: for residuals r ~ N(0, σ²I + λ²11'), Sherman–Morrison and the
  determinant lemma in mpfr supply Σ⁻¹r, Σ⁻¹1, tr Σ⁻¹ and log det Σ to generic
  Gaussian identities: value = -(m log 2π + log det Σ + r'Σ⁻¹r)/2,
  d/dβ = X'Σ⁻¹r, d/dλ = λ((1'Σ⁻¹r)² - 1'Σ⁻¹1),
  d/dσ = σ(|Σ⁻¹r|² - tr Σ⁻¹), mean = λ1'Σ⁻¹r and
  variance = 1 - λ²1'Σ⁻¹1. Checks: 256 against 512 bits (2.5e-67), a dense
  mpfr Cholesky factorization of Σ (8.4e-66), and central differences with
  step 2^-100 at 512 bits for every gradient (2.8e-53).
- Implementation mapping:
  take `case` from the fixture after applying the test helper `fixture_names()`,
  as in `test-linear-marginal.R`; this maps the stored `sd_outcome` name to
  `sd_re_outcome` without changing the frozen fixture.
  `prepared = prepare_data(y ~ x, d ~ w, as.data.frame(case$data), "id", "gaussian")`
  and `out = linear_outcome_marginal(case$theta, prepared, conditional = TRUE)`
  give `out$value`, `out$gradient` (columns in `names(theta)` order, without
  dimnames), `out$mean`, `out$variance`, `out$scale` and the three conditional
  gradient matrices. The selection columns only define the selected rows.
- Coverage: cancellation-free marginal (value, β and σ gradients) and
  the unchanged conditional moments.

## 3. Adaptive inner finite-rule derivative references

The six mixed-derivative constants in `../test-pln-inner.R` were obtained using
mpmath 1.3.0 at 42 decimal digits, with `(a,b,gamma,tau)=(.2,-.3,.8,-.4)`,
`y=4` and five Hermite nodes. They differentiate the moving finite adaptive
rule, not the continuous integral. No package model was called.

For the selected count observation, define
`D=(a+tau*e)/sqrt(1-tau^2)`, `v=b+gamma*e`, and
`h(e)=-e^2/2-log(2*pi)/2+y*v-exp(v)-loggamma(y+1)+log(Phi(D))`, with
`Phi(x)=erfc(-x/sqrt(2))/2`. Generate the unnormalized Hermite rule with
`mp.gauss_quadrature(5, 'hermite')`. For every perturbation, solve `h_e=0`
using `mp.findroot` from `(0,1)`, tolerance `128*mp.eps` and at most 100 steps;
recompute `m` and `H=-h_ee(m)`. Direct multivariate differentiation evaluates

`log(sqrt(2/H) * sum(w_k * exp(h(m + sqrt(2/H)*u_k) + u_k^2)))`.

| Partial | Reference value |
|---|---:|
| `a,a,gamma` | -0.09120888591231401 |
| `a,b,tau` | 0.6758976368282991 |
| `b,b,gamma` | 0.7481695605980433 |
| `a,b` | -0.2155729636950900 |
| `a,b,a` | 0.04542390200999184 |
| `b,a,tau` | 0.6758976368282991 |

The original comparison met scaled tolerance `1e-25`, where the discrepancy
is `abs(computed-reference)/(1+abs(reference))`; its maximum was `1.540e-43`.
This checks distinct coordinates, repeats and a mixed-partial permutation at
these inputs. It does not exhaust all inputs or establish continuous-integral
accuracy. The separate fourth- and fifth-order normal log-CDF constants in
that test were calculated with mpmath 1.3.0 at 65 decimal digits by directly
differentiating `log(erfc(-x/sqrt(2))/2)`.

## 4. Historical saved-object fixtures

These objects were written by historical PanelSelect installations, rather
than reconstructed from current fits. Their `version` metadata records the
source package's version field; the source commit identifies the precise
development snapshot. The neutral schema identifiers describe what the tests
exercise and do not imply separate public releases.

| Fixture | Schema | Source commit | Saved date |
|---|---|---|---|
| `saved-vocabulary.R` | `legacy-control-vocabulary` | `c523ee0597ed39e16b6872a1d2bff44ef16ffaa7` | 2026-10-12 |
| `saved-gk.R` | `gauss-kronrod` | `3bb3b201451190e8c589e0d8cea37f93b16dd0b2` | 2026-10-13 |
| `saved-legacy-settings.R` | `legacy-settings` | `54da535993194fe49d81f5090b43a0da91f84519` | 2026-10-13 |

`saved-vocabulary.R` contains twelve `panel_control()` objects, the
vocabulary-bearing components of six small fits, their data, a summary and a
quadrature check. The fits cover GH Linear selection, staged GL Poisson with
a retained mesh, PLN with a GL observation rule, cross-sectional Linear,
standalone Poisson with a fixed zero SD, and a legacy `probitRE_linearRE()`
fit. `../test-saved-objects.R` rebuilds prepared designs, checks exact preservation
of stored numerical results, and compares independent current refits with a
numerical tolerance. Optimizer effort counts and diagnostic probe choices remain
records of the original run, not cross-platform reproducibility targets. In contrast,
`../helper-vocabulary.R` constructs synthetic historical schemas from current
fits to test normalization; it is not historical-runtime evidence.

`saved-gk.R` contains a Poisson RE fit with GK person effects (8 panels of
3 rows), a PLN RE evaluation with every parameter fixed and GK person and
observation rules, and a `compare = "gk"` check of a Gauss-Hermite fit, with
their data. The check omits per-panel diagnostics because their subnormal
doubles did not round-trip through the original text serialization.
`../test-gk-removal.R` verifies display and the required GL refit guidance.

`saved-legacy-settings.R` contains three controls (default, GL and the checked
binary alias), a Poisson RE fit with Gauss-Hermite person effects, a binary
selection fit with the checked alias, a PLN RE fit from `PLN_RE()` with frozen
adaptive geometry (`fixed_geometry`), and a Linear selection fit from
`probitRE_linearRE()` with a legacy optimizer option (`REPORT`), with their
data. `../test-legacy-removal.R` verifies normalization and refit requirements.

These fixtures use hexadecimal double literals and quoted calls. Formula
environments are restored to `baseenv()` when loaded; the GK and legacy-setting
fixtures also preserve formula attributes explicitly. Their stored values
passed round-trip checks when generated.

## 5. Legacy differential references and preservation

The [PanelSelect local-source reference](panelselect-1.1.2-README.md) and
[PanelCount reference](panelcount-2.0.1.md) document their source identities,
parameter mappings, numerical comparisons and limitations. Source and receipt
hashes are retained in the corresponding fixtures; machine-specific paths are
unnecessary to identify those sources.

All fixtures are frozen. Routine tests require neither historical package
installations nor network access. Never regenerate expected values from the
implementation under test. Any replacement requires independent mathematical
review, an identified reference source, and preservation of the original
numerical evidence.
