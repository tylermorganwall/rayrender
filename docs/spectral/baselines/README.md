# PR 0 Legacy Baselines

These files capture the RGB legacy renderer before the spectral architecture work begins.

- `pr0-legacy-baselines.csv` records seeded render hashes, image sums, render options, CPU time, and peak RSS measurements for representative still and animation scenes.
- `pr0-legacy-environment.md` records compiler, platform, package, and commit information used for the capture.

Regenerate or verify the baseline with:

```sh
tools/spectral-tests/capture-pr0-baselines.R \
  --compare docs/spectral/baselines/pr0-legacy-baselines.csv
```

The image hashes are the stable gate. CPU time and peak RSS are performance baselines and are expected to vary by machine.

