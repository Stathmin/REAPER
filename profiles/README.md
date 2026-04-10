# Snakemake Profiles

RepOrtR ships lightweight Snakemake profile presets:

- `profiles/local/` - local iterative runs.
- `profiles/quiet/` - reduced console noise.
- `profiles/hpc/` - shared cluster defaults.
- `profiles/docker/` - legacy container-oriented defaults.

Usage:

```bash
conda run -n reportr snakemake --profile profiles/local -s Snakefile_modular --use-conda <targets...>
```
