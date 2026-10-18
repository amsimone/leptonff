#!/usr/bin/env python3
"""Build the pinned, optional McMule photonic fixed-order validation client."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from prepare_source import COMMIT, URL, prepare

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, default=ROOT / 'external/mcmule-0.8.0')
    p.add_argument('--build', type=Path, default=ROOT / 'build-mcmule')
    p.add_argument('--jobs', type=int, default=8)
    a = p.parse_args(); source = a.source.resolve(); build = a.build.resolve()
    if not source.exists():
        subprocess.run(['git', 'clone', '--depth', '1', '--branch', 'v0.8.0',
                        '--recurse-submodules', URL, str(source)], check=True)
    patch = prepare(source)
    modules = subprocess.check_output(['git', 'submodule', 'status', '--recursive'], cwd=source, text=True)
    if any(line[:1] != ' ' for line in modules.splitlines()):
        raise ValueError('Missing or unexpected submodule revision')
    subdirs = [line.split()[1] for line in modules.splitlines()]
    trees = {}
    for relative in ['.'] + subdirs:
        directory = source / relative
        if relative != '.' and subprocess.check_output(
                ['git', 'diff', '--name-only', 'HEAD'], cwd=directory, text=True).strip():
            raise ValueError('Modified submodule: ' + relative)
        tracked = subprocess.check_output(['git', 'ls-files'], cwd=directory, text=True).splitlines()
        digest = hashlib.sha256()
        for name in tracked:
            path = directory / name
            if path.is_file():
                digest.update(name.encode() + b'\0' + path.read_bytes())
        trees[relative] = digest.hexdigest()
    env = dict(os.environ, GIT_DIR=str(source / '.git'), GIT_WORK_TREE=str(source))
    meson = shutil.which('meson'); ninja = shutil.which('ninja'); fc = shutil.which('gfortran')
    if not all((meson, ninja, fc)):
        raise RuntimeError('Need meson, ninja and gfortran on PATH')
    build.mkdir(parents=True, exist_ok=True)
    with (build / 'upstream-build.log').open('w') as log:
        if not (build / 'build.ninja').exists():
            subprocess.run([meson, 'setup', str(build), str(source), '-Dgroups=mue',
                            '-Dtests=false', '-Dexamples=false', '--buildtype=release'],
                           env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([ninja, '-C', str(build), '-j' + str(a.jobs)], env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    client_hash = hashlib.sha256((HERE / 'reference.f90').read_bytes() +
                                 (HERE / 'phase_space.f90').read_bytes()).hexdigest()[:12]
    executable = build / ('leptonff_mcmule_' + client_hash)
    rpath = '@executable_path/src' if sys.platform == 'darwin' else '$ORIGIN/src'
    command = [fc, '-O2', '-ffree-line-length-none',
               '-I' + str(build / 'src/libmcmule.dylib.p' if sys.platform == 'darwin' else build / 'src/libmcmule.so.p'),
               '-I' + str(build / 'pre/libmcmule_pre.a.p'), '-I' + str(build / 'src'), '-I' + str(build),
               str(HERE / 'phase_space.f90'), str(HERE / 'reference.f90'), '-J' + str(build),
               '-L' + str(build / 'src'), '-lmcmule',
               '-Wl,-rpath,' + rpath, '-o', str(executable)]
    subprocess.run(command, check=True)
    # Dynamic dependencies are checked before AND after each job. This is not
    # described as a sealed standalone executable: it needs this build tree.
    artifacts = {str(path.relative_to(build)): sha(path) for path in sorted(build.rglob('*'))
                 if path.is_file() and path.suffix in ('.dylib', '.so')}
    artifacts[executable.name] = sha(executable)
    runtime_name = 'libgfortran.dylib' if sys.platform == 'darwin' else 'libgfortran.so'
    runtime = Path(subprocess.check_output([fc, '-print-file-name=' + runtime_name], text=True).strip()).resolve()
    runtime_artifacts = {str(runtime): sha(runtime)} if runtime.is_file() else {}
    manifest = dict(schema='leptonff-mcmule-build-v1', url=URL, tag='v0.8.0', commit=COMMIT,
                    source=str(source), build=str(build), submodules=modules, source_tree_sha256=trees,
                    phase_space_workaround=patch, driver_sha256=sha(HERE / 'reference.f90'),
                    map_sha256=sha(HERE / 'phase_space.f90'),
                    executable=str(executable), artifacts=artifacts, runtime_artifacts=runtime_artifacts,
                    compiler=subprocess.check_output([fc, '--version'], text=True), command=command)
    (build / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(executable)


if __name__ == '__main__':
    main()
