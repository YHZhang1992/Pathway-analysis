# Contributing

1. Create a focused branch and keep changes scoped to one workflow concern.
2. Do not commit credentials, patient-level data, proprietary datasets, or generated run directories.
3. Add or update a synthetic fixture for behavior changes.
4. Run `python -m unittest discover -s tests -v` before opening a pull request.
5. Document changed inputs, outputs, assumptions, and scientific limitations in the README or `docs/`.

Research workflows require numerical and scientific review in addition to passing automated tests.
