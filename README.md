# XRD reference library from the Crystallography Open Database

`build.py` turns COD crystal structures (CIF files) into a powder XRD reference library: for every inorganic entry
it computes the Cu Kα1 pattern with pymatgen and stores the d-spacings, intensities and Miller indices together with
the formula, space group, lattice and a calculated reference intensity ratio (RIR).

The GitHub Action in `.github/workflows/build.yml` downloads COD in 80 slices, builds each slice in parallel and
publishes the merged library as `xrd_library_cod.npz` on the release **xrd-library**.

Left out: COD numbers 9000000 and above (copied from the American Mineralogist Crystal Structure Database), organic
and metal-organic entries (formula with both C and H), and structures with more than 400 atoms per cell.

Data: Crystallography Open Database, https://www.crystallography.net/cod/ (open data). Each phase keeps its COD
number (`id`, e.g. `COD 1010944`) and a link to its COD page.
