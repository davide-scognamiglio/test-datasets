# test-datasets

Data files for [MuSA](https://github.com/davide-scognamiglio/MuSA), fetched by raw URL (MuSA never
clones this repository).

- `dbs_manifest.yaml`: the database manifest `setup` downloads from. Releases pin it by commit SHA.
- `test.vcf.gz`: the single-variant NA12877 VCF used by `-profile test`.
- `get_benchmarking_vcf.sh`: builds a WES-like NA12878 (HG001) callset from the GIAB WGS benchmark VCF.
- `site/`: media for the MuSA website (https://davide-scognamiglio.github.io/MuSA/), kept here so
  the pipeline repository stays small. The website build fetches them by commit SHA and checks
  their SHA-256.
