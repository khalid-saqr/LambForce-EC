# Step 2 Computational Budget

## Execution classes

| Workload | Method | Intended role |
|---|---|---|
| Uniform analytical limits | Lumped model | Unit and regression checks |
| Six-artery primary calculations | Periodic spectral model | Claim-bearing implementation after validation |
| Boundary-condition checks | Bounded finite differences | Verification cases |
| Representative high-fidelity cases | Future 3D model | Later verification only |

## Scaling

For an $`N_x\times N_z`$ periodic grid, the spectral method requires a small number of two-dimensional FFTs and has approximately $`O(N\log N)`$ work with $`N=N_xN_z`$. The finite-difference verification method solves a sparse biharmonic system and is expected to be materially more expensive as the grid grows.

## Frozen budget rule

The publication solver resolution is not inherited from another workflow. Grid size is selected from convergence. Uncertainty calculations use the primary spectral solver; bounded or future 3D models are restricted to predeclared representative verification cases.

## Checkpointing requirement

Step 3 must make artery cases and parameter realizations independently restartable. No scientific result may depend on an uninterrupted notebook session.
