# Public API Audit

This is the intended Capivara 2.0 user-facing surface. The rule is simple:
users see complete science workflows; low-level implementation pieces remain
available to package developers without filling RStudio autocomplete.

## Everyday segmentation

| Function | Why it stays public | Minimal default |
| --- | --- | --- |
| `segment()` | Exact Ward segmentation for manageable cubes. | `Ncomp = 15` |
| `segment_large()` | Sparse-Ward segmentation for real IFU cubes. | `Ncomp = 15`, `knn_k = 40` |
| `build_starlet_mask()` | Explicit starlet support when the user needs it. | Conservative support settings |
| `build_adaptive_support()` | Data-driven alternative support builder. | Adaptive threshold |
| `choose_ncomp_by_snr()` | Decides component count from S/N. | User supplies target S/N |
| `estimate_segment_memory()` | Prevents an accidental exact-Ward memory allocation. | Read-only estimate |
| `summarize_cluster_spectra()` | Region spectra for external fitting. | Flux-aware summaries |
| `reconstruct_cluster_cube()` | Compact representative cube product. | Cluster representatives |
| `reconstruct_flux_preserving_cube()` | Fitting-ready reconstructed cube. | Preserved segment flux |
| `plot_cluster()` / `plot_cluster_spectra()` | The standard visual checks. | Package palette/layout |

## Kinematics

| Function | Scope | Minimal default |
| --- | --- | --- |
| `segment_kinematics()` | Native line maps plus kinematic-aware segmentation. No dynamical assumption. | `knn_k = 50`, `n_segments = 25` |
| `run_kinematic_analysis()` | Kinematic segmentation plus a neutral axisymmetric disc comparison. | `model = "axisymmetric"` |
| `run_manga_bar_model()` | Explicit convenience wrapper for a known barred MaNGA galaxy. | Requires `bar_phi_deg` |
| `kinematic_models()` | Lists installed modelling modules and their requirements. | Read-only registry |

The path-signature segmentation is requested explicitly with
`segmentation_mode = "path_signature"`. Ordinary full-spectrum segmentation is
also explicit through `segment()` or `segment_large()`; it is not silently run
inside a kinematic workflow.

## Expert and Experimental Code

Fitting primitives, raw readers, geometry utilities, diagnostic tables,
batch runners, and component plotters are implementation-level functions. They
are documented for maintainers but should not be the normal package entry
points. The experimental structure scores are also kept out of the headline
API until their scientific interpretation is stable.

## Documentation Policy

Pkgdown remains the website system because it renders R references and
vignettes directly. The site has separate **Kinematic Analysis** and
**Bisymmetric Bar Modelling** articles. Add a future spiral, warp, or outflow
module by registering it in `kinematic_models()`, implementing its model branch,
and adding one dedicated article; this avoids teaching every user a long list of
special-case controls.
