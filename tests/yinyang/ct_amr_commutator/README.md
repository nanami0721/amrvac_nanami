# Yin–Yang CT-AMR curl commutator

The solver initializes the Yin–Yang AMR leaf mesh with a zero magnetic field
and refines a strip of the Yang patch. The test applies the production UCT-HLL
curl operator to scratch face fields for one Faraday update. Targets are
lower-theta boundary faces of both patches over the same local radial and
azimuth intervals, \(r\in[1.18,1.32]\) and \(\phi/(2\pi)\in[0.165,0.225]\).
Each transfer plan verifies that every contributing donor belongs to the
opposite level.

For a receiver face with area \(A\), the diagnostic compares two sequences
for computing its magnetic flux increment \(\Delta\Phi\).

1. Transfer the four native edge EMFs to the receiver and apply the native
   face curl.
2. Apply the native curl on the donor blocks, transfer the resulting face
   field to the receiver, and multiply it by the receiver face area.

Each native edge EMF is prescribed as \(\Delta t\, L E_t\) at the edge
midpoint with \(\Delta t=0.013\). Both increments enter the comparison as
\(\Delta\Phi/(A\Delta t)\). The reference field is the surface average of the
normal component of \(-\nabla\times\boldsymbol E\), evaluated with five-point
Gauss quadrature in each face variable. The field
\(\boldsymbol E=(0.2y^2+0.1zx,-0.3z^2+0.15xy,0.4x^2-0.12yz)\) gives
\(-\nabla\times\boldsymbol E=(-0.48z,0.7x,0.25y)\).

The RMS weights use only the portion of each face inside the fixed radial and
azimuth intervals. The summed area must match the analytic area, including
partially covered faces at the edges of the sample region. The relative RMS
uses the analytic rate RMS over that identical physical support for all three
resolutions. The solver also reports pointwise peak errors and an analytic
comparison for each path. Each rank writes target coordinates, weights, and
commutator rates to its sample file.

The suite doubles all three base-patch cell counts twice, while the level-2
mesh is twice as fine as level 1 in each case. It requires decreasing
commutator RMS at both doublings in each direction, compares the middle
resolution at two and four ranks, checks nontrivial analytic and transfer
signals, and measures the discrete divergence of the receiver curl. The
default geometry uses a wider angular overlap and starts with
\(8\times12\times32\) base cells per patch. The standard geometry starts with
\(8\times16\times48\) cells, using wider azimuth blocks to reduce ghost
storage.

Prepare an external copy of the candidate tree, recording its commit,
working-tree diff, and untracked candidate source files. The runner uses the
project macOS process-tree guard with a 2048 MiB runtime limit.

```sh
python3 tests/yinyang/ct_amr_commutator/run.py EXTERNAL_SOURCE EXTERNAL_RESULTS
python3 tests/yinyang/ct_amr_commutator/run.py EXTERNAL_SOURCE STANDARD_RESULTS --geometry standard
```
