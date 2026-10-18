# Maintainer release checklist

1. Update the CMake project version, `CITATION.cff`, `pyproject.toml`, the Python `__version__`, manual/Doxygen
   release and `CHANGELOG.md` together.
   Document any accuracy, normalisation, ABI or archive compatibility change.
2. Run the standalone CI matrix and check a clean source build. Confirm the
   installed C/C++/Fortran examples and replay executable work after relocation.
3. For changed physics or host code, run the relevant independent references
   and external-generator controls. Public CI does not run external generators
   or certify a differential precision budget.
4. Build the manual with warnings as errors, generate the Doxygen API reference,
   and build/install the Python wheel. Check a CMake build from the Python sdist.
   Review the source distribution and licence notices. Exclude production
   archives, external generator sources, build outputs and credentials.
5. Tag the reviewed commit and publish release notes recording the actual
   validation and remaining limitations. No DOI or precision claim should be
   inferred from a software version tag.

For initial publication of a clean local Git checkout, after authenticating
the GitHub CLI, use:

```sh
gh repo create amsimone/leptonff --public --source . --remote origin --push \
  --description "Lepton fragmentation and matched final-state QED for W/Z observables"
```

Run from the checkout you intend to publish, with a committed `main` branch.
If the target repository already exists, inspect it and configure its remote
instead of creating or replacing it. Repository creation and push do not by
themselves establish that GitHub-hosted CI has passed; check the resulting run.
