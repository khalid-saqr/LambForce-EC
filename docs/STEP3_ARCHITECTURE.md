# LambForce-EC computational architecture

`LambForce-EC` connects the published anisotropic near-wall hydrodynamic field to registered
endothelial membrane mechanics. The published paper remains hydrodynamic ground truth; the package
does not reconstruct Lamb forcing from WSS.

## Hydrodynamic boundary

Version 0.6.0 identifies the actual publication source as the frozen picoNewton v2 notebook and its
scalar/harmonic inputs. Phase 1 independently generates all six arteries in two labelled modes:

- `historical_v2`, for executable regression;
- `verified`, for the mathematically verified real-field nonlinear formulation.

Generated records remain non-claim-bearing until the exact notebook is cold-executed and the all-six
regression is reviewed.

## Mechanics layers

- `published_source.py`: independent published-source hydrodynamic reproduction;
- `models.py`: artery, geometry, material, and run contracts;
- `loads.py`: signed integration, isotropic/anisotropy controls, and force-conserving distributions;
- `solvers/spectral.py`: periodic Fourier-spectral primary mechanics solver;
- `solvers/bounded_fd.py`: bounded verification solver;
- `solvers/lumped.py`: analytical foundation baseline;
- `workflow.py`: simultaneous normal Lamb loading and tangential WSS mechanics;
- `protocol.py`: source, parameter, traceability, and claim-bearing gates;
- `io.py`: checksummed standardized inputs and result manifests.

No fitted transfer efficiency, arbitrary equivalent-pressure gain, or outcome-selected localization is
permitted.
