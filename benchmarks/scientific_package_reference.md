# Scientific package refactor validation (0.26)

The refactor separates a standalone fragmentation/evolution library from the
full matched W/Z package, retaining the original full archive for generator
relinkers. Special functions and the x-space decay validation helper were
moved without changing their formulas. Existing kernels, threshold maps,
matching prescriptions, signed production records and covariance are retained.

Local validation on macOS ARM:

| Configuration | Checks | Outcome |
|---|---:|---|
| AppleClang shared, Fortran module enabled, all research references | 77 | Passed |
| GCC static, Fortran module enabled, research executables disabled | 42 | Passed |
| Standalone static core | 11 | Passed |
| Shared core built from the Python source distribution | 12 | Passed |

The shared suite includes the previous analytic/independent numerical checks,
W/Z matching and production replay regression tests, the new core C API,
Python/native lifecycle/error checks and the canonical analytic fixture.
The installed-consumer check builds C++, C and Fortran programs using only a
relocated installation and installed examples. Core-only consumers additionally
verify that requesting an unavailable matched component fails. Final packaging
smoke checks cover runtime version reporting, installed reference data and
pkg-config linking where the tool is available.

The HTML manual builds with Sphinx/MyST warnings treated as errors. The Python
wheel and source distribution build successfully; the installed wheel runs a
native fragmentation calculation. The source archive contains the native
sources, examples, tests, canonical data and licence notices, while production
archives and external generator installations remain excluded.

CI covers GCC static/shared and Clang shared full suites, a standalone shared
core, and manual/Python-distribution builds. The documentation/distribution job
also checks a native core built from the sdist and generates the Doxygen API
reference. Local and hosted checks are separate evidence; consult the relevant
GitHub Actions run for hosted status.

These checks establish package/interface behavior and preserve existing
numerical results. They do not certify new per-bin phenomenological precision,
provide additional high-statistics QCD production, implement native DYTurbo
matching, or supply portable interpolated evolution/decay-response tables.
