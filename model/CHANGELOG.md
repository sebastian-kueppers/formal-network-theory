# Changelog for `cooccurrence.model`

All notable changes to this function are documented here.  
This file follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) principles.

## [1.1.3] - 2026-11-10
### Changed
- **Performance:** the Glauber dynamics are vectorised over all goal activations. Instead of looping over activations and, within each, over `n.Glauber` iterations, all activations run as independent chains (rows of a K × N state matrix) that are updated simultaneously. In each iteration, every chain updates one randomly chosen node (random-site Glauber dynamics), as before.
	- The flip probability is computed from the local energy change of the chosen node, ΔH = 2·s_n·(τ_n + Σⱼ W_nj s_j), instead of evaluating the full Hamiltonian twice. This is algebraically identical (zero diagonal) and reduces the cost per update from O(N²) to O(N).
	- The full chain (`n.Glauber` × N per activation) is no longer stored. Node states are summed only within each activation's averaging window.
	- The 50-iteration averaging window is drawn per activation before the dynamics start, from the same range as before (iterations `n.ignored.Glauber + 1` to `n.Glauber - 49`). Iterations after the last window are not computed.
- **Results are equivalent in distribution to 1.1.2, but not bit-identical for the same seed**, because random numbers are drawn in a different order. Simulations from ≤ 1.1.2 cannot be reproduced exactly with 1.1.3.
- The progress bar now counts Glauber iterations instead of goal activations.
- `stateE1`/`stateE2` are computed with `rowSums(..., drop = FALSE)` instead of `apply(..., 1, sum)`. The result is the same, but it no longer fails if an emotion has only one node.

### Removed
- Internal function `hamiltonian()`, which is no longer needed.

### Unchanged
- Arguments, defaults, input checks, and the structure and column names of the returned list (`weights.matrix`, `initial.thresholds`, `data`, `data.thresholds`, `network.coding`).

## [1.1.2] - 2026-10-12
### Added
- Fixed thresholds for existing networks. If `threshold.E1` and `threshold.E2` are both supplied with `emotion.model = "existing.network"`, they replace the network's thresholds, and none are looked up in or converted from the network.
	- Values are interpreted on the {-1,1} scale used by the Glauber dynamics and are not converted. Only the weights are still divided by 4 if the network was estimated on {0,1} data, so `network.coding` is still required for plain matrices.
	- Nodes only in `names.E1` get `threshold.E1`, nodes only in `names.E2` get `threshold.E2`, and shared nodes get their mean (same rule as for generated networks).
	- Supplying only one of `threshold.E1`/`threshold.E2` with an existing network stops with an error.
	- Combining fixed thresholds with `network.thresholds` stops with an error.
	- The out-of-range warning (|threshold| > 1) and the range checks on `threshold.E1`/`threshold.E2` also apply to fixed values.
- `ising_network()`: `thresholds` is now optional (default `NULL`), so `ising_network(weights, coding = ...)` can be used together with fixed thresholds. Without fixed thresholds, a network lacking thresholds still stops with the "No thresholds found" error.

### Changed
- **Behavior change (supersedes 1.1.0):** `threshold.E1`/`threshold.E2` are no longer always ignored for existing networks. They are ignored only if not both are explicitly supplied; then the network's thresholds are used as before.
- Docs for `threshold.E1`, `threshold.E2`, `network.thresholds`, `network.coding` and `ising_network()` updated accordingly.

### Removed
- The warning "threshold.E1 and threshold.E2 are ignored" for existing networks, since it would fire exactly in the new fixed-threshold use case.

### Fixed
- Integer check for `n.nodes.in.network.E1`/`.E2` used `||` inside `!()`, so it only failed if both were non-integer. It now fails if either is (`&&`). Typo in the error message corrected.
- `dist.goal.activation` sum-to-1 check compared floating-point numbers with `!=` and could spuriously fail (e.g., for `c(0.3, 0.3, 0.4)`). It now uses `all.equal()`.
- Roxygen defaults for `n.shared.nodes` (3) and `n.cross.neg.edges` (8) now match the function signature.

## [1.1.1] - 2026-10-07
### Added
- New exported helper `ising_network(weights, thresholds, coding)`. The helper bundles edge weights, node thresholds and the coding of the estimation data (`"0/1"` or `"-1/1"`) into an object of class `"ising_network"`, independent of the estimation software (glmnet, IsingFit, IsingSampler, bgms).
	- `coding` has no default and must be given explicitly. This prevents silent mis-scaling of {0,1} estimates.
	- Validates numeric, symmetric weights with identical row/column names and a zero diagonal. Thresholds are numeric without `NA`. They are matched to nodes by name, never by position, or taken in column order if unnamed.
- `existing.network` now accepts an `ising_network` or a plain list with the elements `weights`, `thresholds` and (optionally) `coding`. A list without `weights` or `thresholds` stops with an informative error.

### Changed
- `network.coding = "auto"` now takes the coding from the `coding` element of an `ising_network`. If `network.coding` is set explicitly and contradicts the stored coding, the function stops.
- Error messages for missing thresholds and uninferrable coding now point to `ising_network()`.
- The out-of-range threshold warning now says "suggested" instead of "allowed". It is a warning, not a constraint.
- Docs for `existing.network` and `network.coding` updated. Matrices and qgraph objects are documented as requiring `network.thresholds` and `network.coding`.

### Removed
- The `message()` summarizing max |weight|, threshold range and input coding for `existing.network`. The out-of-range warnings remain.

---

## [1.1.0] - 2026-10-06
### Changed
- `emotion.model = "existing.network"` now requires an Ising network, because the Glauber dynamics need thresholds as well as weights.
	- bootnet results estimated with non-Ising defaults (EBICglasso, pcor, cor, huge, ggmModSelect, TMFG, LoGo, GGMncv, adalasso, relimp) are rejected.
	- Thresholds come from the network (`bootnetResult$intercepts`, or attribute `"intercepts"` of a matrix) or from `network.thresholds`.
	- **Behavior change:** `threshold.E1`/`threshold.E2` are ignored for existing networks, with a warning if supplied. The thresholds of the supplied network are used instead.
- Weights are validated: numeric, named, symmetric, zero diagonal.

### Added
- `network.thresholds`: numeric vector (named by node, or in column order of the network) or one-column matrix.
- `network.coding` (`"auto"`, `"0/1"`, `"-1/1"`). `"auto"` infers the coding from `bootnetResult$data`. Weights and thresholds fitted on {0,1} data are converted to {-1,1} (Haslbeck, Epskamp, Marsman & Waldorp, 2021, Table 2): β = β*/4, α = α*/2 + Σⱼβ*ᵢⱼ/4.
- Warnings if converted weights exceed |0.8| or thresholds exceed |1|.
- A message summarizing max |weight|, threshold range and input coding.

---

## [1.0.3] - 2026-08-18
### Fixed
- Feasibility checks (shared-node count, central cross-edge count, random.network edge maximum, small-world wiring) are now skipped when emotion.model = "existing.network". These parameters don't apply to imported networks, so validating them previously produced irrelevant errors/warnings.
- Spurious "weight parameters are ignored" warning under emotion.model = "existing.network". Internal reassignment of weight.pos.E1.E2/weight.neg.E1.E2 made missing() report them as user-supplied, firing the warning on every call. Supplied-argument status is now read from match.call() before any reassignment.

---

## [1.0.2] - 2026-08-13
### Fixed
- Shared nodes could not exceed each network's unique-node count. Networks were generated at unique size (`n.nodes.unique.E1`/`.E2`) and shared-node connections were then grafted on via a per-shared-node loop that reused unique-node rows as templates (`mask1[i, ]`). This produced an out-of-bounds error (`subscript out of bounds`) whenever `n.shared.nodes > n.nodes.unique`, and even below that threshold assigned shared-node edges from an arbitrary node correspondence rather than the generated topology.
- Simulated single-emotion activation (`kind == "E1"`/`"E2"`) drove shared nodes twice — once as members of `idx.E1`/`idx.E2` (which already include the shared nodes) and again via a separate `idx.shared` step — so shared nodes received `+2·intensity` while all other nodes received `+intensity`. Removed the redundant `idx.shared` activation; shared nodes are now driven once, like every other node. This bias was negligible at the small `n.shared.nodes` reachable in ≤1.0.1 and becomes consequential only at the high shared-node fractions 1.0.2 enables.

### Changed
- Networks are now generated at their **full** size (`n.nodes.in.network.E1`/`.E2`, including shared nodes) across all generated model types (`small.world.network`, `preferential.attachment.network`, `random.network`, `fully.connected.network`). Shared nodes are therefore real nodes carrying the generated topology, rather than being reconstructed after the fact.
- The combined weight matrix is now assembled by direct block placement of each network's sub-matrices (unique–unique, unique–shared, shared–unique) instead of the previous `rbind`/`mask` loop. This removes the node-by-node loop entirely and preserves each model's topology for any `n.shared.nodes` up to `min(n.nodes.in.network.E1, n.nodes.in.network.E2) - 1`.
- Shared–shared edges now follow the union of E1's and E2's generated topologies, with the edge weight averaged where both networks connect the pair. (Previously the shared block was set to a fully connected `(weight.E1 + weight.E2) / 2`.)

### Added
- Input validation for the new full-size regime:
	- `stop` if `n.shared.nodes >= n.nodes.in.network.E1` or `>= n.nodes.in.network.E2` (each network must retain at least one unique node).
	- `stop` if `cross.edge.type = "central"` and the number of cross edges exceeds `min(n.nodes.unique.E1, n.nodes.unique.E2)` (previously produced silent `NA` indices).
	- `stop` if requested `random.network` edges exceed the complete-graph maximum `n*(n-1)/2`.
	- `warning` if `n.nodes.in.network` is small relative to `sw.nei`, where the small-world ring may be degenerate.

---

## [1.0.1] - 2026-07-29
### Changed
- replaced `abs.weight.E1.E2` (a single magnitude applied to all cross edges) with two separate parameters: `weight.pos.E1.E2` for positive cross edges (a positive value, default 0.08) and `weight.neg.E1.E2` for negative cross edges (a negative value, default -0.08).
	- wrong-signed inputs are coerced to the correct sign with a warning (e.g., a positive `weight.neg.E1.E2` is changed to its negative).
	- the `> 0.8` magnitude cap and the `existing.network` "ignored argument" warning were updated to cover both new parameters.

### Fixed
- `random.network` generated network 2 with E1's unique-node count (`n.nodes.unique.E1`) instead of E2's, causing a matrix-assembly error whenever emotion 1 and 2 had different sizes. Now uses `n.nodes.unique.E2`.

---

## [1.0.0] - 2025-11-28
### Added
- User can choose different network model types for network creation: `small.world.network`, `preferential.attachment.network`, `random.network`, `fully.connected.network`, `existing.network`
- Support for passing an existing network via `existing.network`.
	- Node names for networks `names.E1` and `names.E2` must be provided when using an existing network.
- New parameters for network structure
	- `r.m.components.times.x` (for random networks).
	- `sw.nei`, `sw.p` (small-world network parameters).
	- `pa.power`, `pa.m` (preferential attachment network parameters).
- Cross edges can be chosen to be random or attached to the most central node in each network, with `cross.edge.type`.
- Roxygen2 documentation added for all function arguments and return values.

### Changed
- renamed `n.positive.cross.edges` / `n.negative.cross.edges` to `n.cross.pos.edges` / `n.cross.neg.edges`.
- boolean `shared.node` removed; instead the user can set any number of shared nodes, `n.shared.nodes` 