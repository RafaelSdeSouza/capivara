# Reproducing the Paper

The published paper is **de Souza et al. (2025), “capivara: a
spectral-based segmentation method for IFU data cubes”, MNRAS 539,
3166–3179**, DOI
[10.1093/mnras/staf688](https://doi.org/10.1093/mnras/staf688), arXiv
[2410.21962](https://arxiv.org/abs/2410.21962).

## Published vs ongoing work

The 2025 article is the canonical reference for the published Capivara
method. The repository also contains ongoing Capivara 2.0 research,
tutorial scripts, and native kinematic interfaces. Those later materials
must not be described as part of the published 2025 result unless a
file-level provenance mapping demonstrates that relationship.

## Reproduction status

The repository includes existing real-data panels and research scripts,
but a complete mapping from each published figure to an input cube,
access identifier, configuration, random seed, generator, and
machine-readable result table has not yet been verified. Consequently,
this guide does not offer a single command that claims to reproduce the
paper.

The website provenance table records what is known and leaves unresolved
fields explicit. A complete reproduction entry should identify:

1.  the MaNGA input cube and access identifier;
2.  wavelength sampling and spatial/WCS handling;
3.  support construction;
4.  representation and clustering backend;
5.  all parameter values and random seeds where applicable;
6.  the executable generator;
7.  the machine-readable table or saved result object; and
8.  the final figure path.

## Citation

Use the package citation as the single machine-readable record:

``` r
citation("capivara")
```

The citation includes all ten published authors, the journal, volume,
issue, page range, DOI, and arXiv identifier.
