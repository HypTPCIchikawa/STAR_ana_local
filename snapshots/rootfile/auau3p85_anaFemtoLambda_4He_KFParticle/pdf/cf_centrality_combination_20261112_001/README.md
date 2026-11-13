# Centrality combination comparison

This diagnostic uses the verified saved Lambda–4He histogram export, selecting only six selected-signal QA histograms. It sums all nSigma y bins including flow and keeps regular k* bins. It does not reopen the production ROOT or change event pairing.

`global` reproduces the existing CF: sum raw SE and ME and normalize once. `balanced2` and `balanced3` match each partition ME to its SE yield in the diagnostic 0.4–0.5 GeV/c interval before combining counts. They are inclusive ratios, equivalently averages of individual CFs with k-dependent scaled-ME weights. `fixed_norm_SE_weights2/3` instead average individual CFs with constant SE normalization-yield fractions. These are different observables; they are not expected to agree in general. No physical plateau is asserted for the matching interval.

The five methods are stored for high-k window-integral ratios and 0.020 GeV/c displayed bins. Window ratios are not arithmetic averages of CF points. The uncertainty uses stored cell variances and analytic first derivatives including estimated normalizations and estimated fixed weights. Differences use the shared input covariance. Event reuse, mixed-pair correlations and other inter-cell covariance are unavailable.

Run `python3 source/compare_combination.py settings.json` in a fresh copy of this study. The runner refuses existing results. Configuration, source, raw projected counts, normalization counts, window results, displayed points and validation with fingerprints are retained. No images are generated.

Validation includes the original-export SHA, unchanged input fingerprint, previous global-window values, normalization-domain closure, equality of the three inclusive methods when all class scales are equal, and a toy where both individual CFs are one but unequal ME exposures distort an unbalanced raw sum. The fixed-weight mean is explicitly excluded from the common-scale identity because it represents a different estimand.

Normalization closure for the fixed-weight estimator means treating H as one window: every individual C_i(H) is one, so their average is one. It does not assert that a D(k)-weighted integral of the averaged point sequence is exactly one. The latter is recorded separately without rescaling; no additional global normalization is applied to fixed-weight averages.
