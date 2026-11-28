# SeekDeep documentation site

Source for <https://bailey-lab.github.io/SeekDeep/>, built with [Quarto](https://quarto.org/).
`.github/workflows/docs.yml` builds it on pull requests that touch `docs/`, and builds and
publishes it on pushes to `master` and on `v*` tags.

## Building locally

```bash
quarto preview docs      # live preview while editing
quarto render docs       # full build into docs/_site/
```

No R or Python is needed to render: pages are plain Markdown.

## Command reference

`reference/` is generated from a SeekDeep build. Don't edit it by hand. After changing flags or
commands, rebuild SeekDeep, regenerate the reference and commit the result:

```bash
python3 docs/_scripts/generate_reference.py                     # uses ./bin/SeekDeep
python3 docs/_scripts/generate_reference.py --seekdeep /path/to/SeekDeep
```

Regenerate before tagging a release, so the published reference matches the tagged code.

## Keeping the site small

The CI build fails if any file is over 2 MB or the site is over 30 MB. Put example data and other
large downloads on Zenodo or in release assets and link to them. Export images at web size
(about 1600-2400 px wide).
