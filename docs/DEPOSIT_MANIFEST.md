# Deposit manifest

Concept DOI: 10.5281/zenodo.22693563 (always resolves to the latest version)

## Repo path to deposited filename

| repo path | deposited filename | first deposited in |
|---|---|---|
| `docs/composition_note.md` | `SSc_composition_note.md` | v1, 10.5281/zenodo.22693564 |
| built from `docs/composition_note.md` | `SSc_composition_note.pdf` | v1, 10.5281/zenodo.22693564 |
| `figures/composition_fig.png` | `composition_fig.png` | v1, 10.5281/zenodo.22693564 |

Filenames differ between the repo and the deposit for historical reasons. They are
mapped here rather than renamed, so that version diffs stay readable.

## Figure

`src/make_composition_fig.py` generates both `figures/composition_fig.png` and
`figures/composition_fig.pdf`. **The png is the deposited form**, because the note is
written in markdown and references it by name. The pdf is the vector form, kept in the
repo for reuse at print size and not deposited.

The deposited figure is therefore reproducible from this repo, which is what the note's
Data and code section claims.

## Build

The PDF is generated from the markdown and is never edited directly:

```
pandoc docs/composition_note.md -o SSc_composition_note.pdf \
  --pdf-engine=xelatex -V geometry:margin=1in -V fontsize=11pt
```

The markdown references `composition_fig.png` by bare filename, so the png must be
resolvable from the working directory at build time, or passed via `--resource-path`.

## Versions

| version | DOI | change |
|---|---|---|
| v1 | 10.5281/zenodo.22693564 | initial deposit, 2026-09-10 |
| v2 | (on publication) | corrected the first author of the reanalyzed study (Jarnagin, not Ashton) and its publication year (2026, not 2025). No data, method, figure or result changed. |
