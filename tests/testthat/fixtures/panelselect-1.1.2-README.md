# PanelSelect 1.1.2 differential fixture

`panelselect-1.1.2.R` contains outputs executed by a preserved local
PanelSelect 1.1.2 source snapshot, frozen on 2026-09-26. This is not asserted
to be an unmodified CRAN archive. The fixture records its source hashes, Git revision, rebuilt archive
hash, runtime receipt hashes, inputs and parameter mappings. Runtime tests need
neither an older PanelSelect installation nor PanelCount.

Sixteen regular cases cover binary/Gaussian selection, two correlation signs,
all four rho/tau free/fixed-zero combinations, unequal nonunit scales, irregular
panels, interleaved input and all-unselected panels. H8 rules share their exact
finite objective: the general Gaussian route when tau is free and the common
conditional outcome reduction when tau is fixed zero. Reviewed development
comparisons also retained H16 and H32, two finite-difference steps and the
alternative Gaussian grid. These larger runs are excluded from routine tests.

Natural per-panel scores are mapped from the legacy log-SD/Fisher coordinates.
Prediction references come from executed original helpers with the same named
natural covariance: selection probability, potential outcome, their gradients
and row SEs, and the potential-outcome average and its gradient-based SE.
The public modern prediction tests omit response columns.

One deliberate tail case exposes the old 1e-100 probability floor. Its old
likelihood and zero gradient are retained as evidence of the discrepancy;
the modern result must instead match the analytic independent-normal/probit
convolution and its nonzero gradient.

The fixture retains relative source names and their SHA256 values, Git revision
`9a15cc8e2dea68810b440cef5f20775a9505499c`, comparison-receipt identifiers and
hashes, and rebuilt source archive SHA256
`1cb3bc9f38c0d206ca3af29969bcdb0161f8cc4eaf6855eabdd1e60f1a6e45ee`.
These identify the original evidence without depending on a local path.

The text fixture preserves the original saved numerical values. Regeneration
would require the identified source and an isolated historical installation,
independent review of mappings and outputs, and a separate replacement fixture.
Routine tests do not regenerate these values. Never refresh the expected values
from the modern evaluator.
