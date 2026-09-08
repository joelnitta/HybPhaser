# Changelog

## rhybphaser 0.1.0.9000

- Renamed the package from **HybPhaser** to **rhybphaser** to avoid
  confusion with the original script collection at
  <https://github.com/LarsNauheimer/HybPhaser>. This is an R-package
  reimplementation of the HybPhaser method by Nauheimer et al.
- Unified how external tools are executed.
  [`run_generate_consensus_sequences()`](https://joelnitta.github.io/rhybphaser/reference/run_generate_consensus_sequences.md),
  [`run_extract_mapped_reads()`](https://joelnitta.github.io/rhybphaser/reference/run_extract_mapped_reads.md),
  [`run_clade_association()`](https://joelnitta.github.io/rhybphaser/reference/run_clade_association.md)
  and
  [`run_phasing()`](https://joelnitta.github.io/rhybphaser/reference/run_phasing.md)
  (and their `*_from_config()` variants) now take
  `engine = c("docker", "local")`, defaulting to `"docker"` (the pinned
  `rhybphaser` image). This replaces the previous inconsistency where
  consensus generation was Docker-only while BBSplit ran locally first
  with a Docker fallback. The `docker_fallback` argument is removed; use
  `engine = "local"` to run tools from `PATH`.
- The Docker image was renamed `joelnitta/hybphaser` -\>
  `joelnitta/rhybphaser` and is published from CI
  (`.github/workflows/docker.yml`).
- [`hybpiper_stats()`](https://joelnitta.github.io/rhybphaser/reference/hybpiper_stats.md)
  now falls back to a native R implementation when `hybpiper stats`
  fails in `"supercontig"` mode (works around a `KeyError` in HybPiper
  2.3.4 with gene + OG style target headers). The fallback writes
  HybPiper-compatible `seq_lengths.tsv` and `hybpiper_stats.tsv` from
  the assembled `*_supercontig.fasta` files.
