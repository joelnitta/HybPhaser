# rhybphaser 0.1.0.9000

* Renamed the package from **HybPhaser** to **rhybphaser** to avoid confusion
  with the original script collection at
  <https://github.com/LarsNauheimer/HybPhaser>. This is an R-package
  reimplementation of the HybPhaser method by Nauheimer et al.
* `hybpiper_stats()` now falls back to a native R implementation when
  `hybpiper stats` fails in `"supercontig"` mode (works around a `KeyError` in
  HybPiper 2.3.4 with gene + OG style target headers). The fallback writes
  HybPiper-compatible `seq_lengths.tsv` and `hybpiper_stats.tsv` from the
  assembled `*_supercontig.fasta` files.
