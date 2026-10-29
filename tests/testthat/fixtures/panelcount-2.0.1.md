# Executed PanelCount 2.0.1 numerical reference

`panelcount-2.0.1.R` contains actual outputs of the unchanged original package,
compiled and installed into a separate library from PanelSelect 2.0.0. Its 27 input
file hashes are embedded. All 17 R/native/metadata files in the preserved baseline
matched the original byte-for-byte. The fixture records the source archive,
installed-runtime receipt, recipe and worker hashes, R session, and FD step sizes.
The original source tarball SHA256 is
`7662f8c35c11910ae0faccd3733e0f4ff37ef64c17e357b023aa62d17669979d`.

Eight cases cover the five models: standalone probit, Poisson, PLN (both its
ordinary fixed and explicitly frozen-prior code paths), and the two selected
count models at zero and nonzero correlations. Each has four three-row panels,
including an all-unselected panel for selection models. Natural parameters,
nonunit scales, data, supported offsets, and Hermite orders are explicit.

Matching-rule likelihood comparisons at order 12 are exact up to floating-point
arithmetic. Scores are compared to finite differences of the actual old
objective at steps 5e-5 and 1e-5, excluding free correlation columns whose tensor
grids move in different orientations. Raw old scores are retained: standalone
Poisson/PLN multiply SD derivatives by an extra scale factor, so those raw scores
are not correctness oracles. At nonzero correlations, the different old/new
tensor orientations require order-refined approximate comparisons; the saved
order 30 cases use tolerances justified by a separately retained 12/20/30 ladder.
These comparisons do not claim equal finite-grid correlation derivatives or
certify general integration accuracy.

Prediction references are actual calls to the original two public helpers.
`prob` maps to selection probability, `outcome` to the potential mean, and
`actual_outcome` to the joint mean. Modern public `predict()` runs on explicit
fixed-parameter objects, and both modern compatibility prediction helpers run on
the corresponding legacy parameter vectors. No standalone public prediction helper exists in the
original package; no synthetic substitute is labeled as executed legacy code.

The fixture retains relative source filenames and hashes, the source archive
hash above, runtime-receipt, worker and recipe hashes, R session information,
and finite-difference step sizes. These identify the executed reference without
depending on private file locations.

The text fixture preserves the original saved numerical values. Any replacement
requires the identified original source, an isolated historical installation,
independent review of mappings and outputs, and preservation of the original
reference. Routine tests read the fixture and call modern code; they do not
install PanelCount, run a legacy optimizer, require network access, or regenerate
reference values. Never refresh expectations from the modern evaluator.
