#!/usr/bin/env python3
"""Build pinned DYTurbo 1.4.2 and its native client in isolated directories."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import time

from prepare_source import prepare

ROOT = Path(__file__).resolve().parents[2]
# Digests of the official archives actually used for the native validation.
ARCHIVES = {
    'dyturbo-1.4.2.tar.gz': ('https://dyturbo.hepforge.org/downloads/dyturbo-1.4.2.tar.gz',
        '91d4ab97966a5dc88d7d63a0588bbc16a946980c57b35ae802c62d8935f94b98'),
    'LHAPDF-6.5.4.tar.gz': ('https://lhapdf.hepforge.org/downloads/LHAPDF-6.5.4.tar.gz',
        '2443a4b32cc3b0597c8248bd6e25703ace9c91a7a253c5f60b1b5428ef9c869e'),
    'Cuba-4.2.2.tar.gz': ('https://feynarts.de/cuba/Cuba-4.2.2.tar.gz',
        '8d9f532fd2b9561da2272c156ef7be5f3960953e4519c638759f1b52fe03ed52'),
    'cubature-1.0.3.tar.gz': ('https://codeload.github.com/stevengj/cubature/tar.gz/refs/tags/v1.0.3',
        'a4788fc67dc66639839e73738dc2c1df81f3a3461d35e99b4a4a49b3a9deda10'),
    'HELLx.v3.0.tgz': ('https://www.roma1.infn.it/~bonvini/hell/downloads/HELLx.v3.0.tgz',
        '0c424f09091c904683166a74e98cd1c1e53e2d83eea5fbfb16d65181b386642e'),
    'HELLx-data.v3.tgz': ('https://www.roma1.infn.it/~bonvini/hell/downloads/HELLx-data.v3.tgz',
        'f456f9d9bc34dc3b01dcac6c0cb020cd8502294e594daa4549e291b2f65e97ca'),
}


def digest(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def fetch(cache, name):
    url, expected = ARCHIVES[name]
    path = cache / name
    if not path.exists():
        partial = cache / (name + '.partial')
        subprocess.run(['curl', '-L', '--fail', '--max-time', '180', '--output', str(partial), url], check=True)
        if digest(partial) != expected:
            raise RuntimeError('Unexpected archive content: ' + name)
        partial.rename(path)
    if digest(path) != expected:
        raise RuntimeError('Archive checksum mismatch: ' + name)
    return path


def extract(archive, destination, root):
    """Extract regular files; materialise safe internal archive links as copies."""
    destination = destination.resolve()
    with tarfile.open(archive, 'r:gz') as t:
        members = t.getmembers()
        if len(members) > 50000 or sum(x.size for x in members) > 2 * 1024**3:
            raise RuntimeError('Unexpected dependency archive size')
        for m in members:
            p = Path(m.name)
            target = destination / p
            if (p.is_absolute() or '..' in p.parts or not p.parts or p.parts[0] != root
                    or not target.resolve().is_relative_to(destination)
                    or not (m.isdir() or m.isfile() or m.issym() or m.islnk())):
                raise RuntimeError('Unsafe dependency member: ' + m.name)
            if m.issym() or m.islnk():
                link = Path(m.linkname)
                linked = (target.parent / link) if m.issym() else (destination / link)
                if link.is_absolute() or not linked.resolve().is_relative_to(destination / root):
                    raise RuntimeError('Escaping dependency link: ' + m.name)
        for m in members:
            target = destination / m.name
            if m.isdir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with t.extractfile(m) as src, target.open('wb') as dst:
                    shutil.copyfileobj(src, dst)
                target.chmod(m.mode & 0o777)
                # Generated autotools files must retain the distribution
                # timestamps; extraction order must not trigger autoreconf.
                os.utime(target, (m.mtime, m.mtime))


def replace_once(path, before, after):
    data = path.read_text()
    if after in data:
        return
    if data.count(before) != 1:
        raise RuntimeError('Unrecognised upstream source at ' + str(path))
    path.write_text(data.replace(before, after))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work', type=Path, default=ROOT / 'external/dyturbo-native')
    p.add_argument('--cache', type=Path, default=ROOT / 'external/downloads')
    p.add_argument('--client-build', type=Path, default=ROOT / 'build-dyturbo-interface')
    p.add_argument('--cc', default='clang' if os.uname().sysname == 'Darwin' else 'cc')
    p.add_argument('--cxx', default='c++')
    p.add_argument('--fc', default='gfortran')
    p.add_argument('--jobs', type=int, default=4)
    a = p.parse_args()
    if a.jobs < 1:
        p.error('--jobs must be positive')
    work, cache = a.work.resolve(), a.cache.resolve()
    work.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    archives = {name: fetch(cache, name) for name in ARCHIVES}
    logdir = work / 'logs'
    logdir.mkdir(exist_ok=True)
    commands = []
    env = dict(os.environ)

    def run(label, args, cwd):
        print(label, flush=True)
        commands.append(dict(stage=label, argv=[str(x) for x in args], cwd=str(cwd)))
        with (logdir / (label + '.log')).open('w') as log:
            subprocess.run([str(x) for x in args], cwd=cwd, env=env, check=True,
                           stdout=log, stderr=subprocess.STDOUT)

    lha = work / 'LHAPDF-6.5.4'
    if not lha.exists():
        extract(archives['LHAPDF-6.5.4.tar.gz'], work, lha.name)
    run('lhapdf-configure', ['./configure', '--prefix=' + str(lha / 'local'),
        '--disable-python', '--enable-lhaglue', 'CC=' + a.cc, 'CXX=' + a.cxx], lha)
    run('lhapdf-build', ['make', '-j' + str(a.jobs)], lha)
    run('lhapdf-install', ['make', 'install'], lha)
    env['PATH'] = str(lha / 'local/bin') + os.pathsep + env['PATH']

    prepared = work / 'source'
    source = prepared / 'dyturbo-1.4.2'
    if not prepared.exists():
        prepare(archives['dyturbo-1.4.2.tar.gz'], prepared)
    manifest = json.loads((prepared / 'source_manifest.json').read_text())
    if manifest['archive_sha256'] != ARCHIVES['dyturbo-1.4.2.tar.gz'][1]:
        raise RuntimeError('Prepared source belongs to a different archive')
    cuba = source / 'Cuba-4.2.2'
    if not cuba.exists():
        extract(archives['Cuba-4.2.2.tar.gz'], source, cuba.name)
    # Cuba's archive-member make rules race under parallel make on macOS.
    run('cuba-configure', ['./configure', 'CC=' + a.cc, 'CFLAGS=-O2 -fPIC -fcommon'], cuba)
    run('cuba-build', ['make', '-j1'], cuba)
    with tarfile.open(archives['cubature-1.0.3.tar.gz'], 'r:gz') as t:
        for name in ['hcubature.c', 'pcubature.c', 'cubature.h', 'clencurt.h', 'converged.h', 'vwrapper.h']:
            (source / 'cubature' / name).write_bytes(t.extractfile('cubature-1.0.3/' + name).read())
    if not (source / 'HELLx/Makefile').exists():
        extract(archives['HELLx.v3.0.tgz'], source, 'HELLx')
    if not any((source / 'HELLx/data/LL').glob('*.table')):
        extract(archives['HELLx-data.v3.tgz'], source, 'HELLx')
    matrix = source / 'HELLx/include/math/matrix.hh'
    replace_once(matrix, 'det(*this);', 'det();')
    replace_once(matrix, 'd==0)', 'd==0.)')
    run('hellx-build', ['make', '-j' + str(a.jobs), 'lib', 'CC=' + a.cxx,
                      'CFLAGS=-O2 -fPIC -Wall -I./include'], source / 'HELLx')
    # Equivalent zero initialisation; Clang rejects initialised VLAs in C++.
    ct = source / 'counterterm/ctint.C'
    for name in ['sig11_qed', 'sig12_qed']:
        replace_once(ct, 'double ' + name + '[mesq::totpch] = {0.};',
                     'double ' + name + '[mesq::totpch];\n      for (int channel=0; '
                     'channel<mesq::totpch; ++channel) ' + name + '[channel]=0.;')
    libname = 'libgfortran.dylib' if os.uname().sysname == 'Darwin' else 'libgfortran.so'
    lib = Path(subprocess.check_output([a.fc, '-print-file-name=' + libname], text=True).strip()).resolve()
    if not lib.is_file():
        raise RuntimeError('Cannot locate the Fortran runtime: ' + str(lib))
    ldflags = '-L' + str(lib.parent) + ' -Wl,-rpath,' + str(lib.parent)
    prefix = source / 'local'
    run('dyturbo-configure', ['./configure', '--prefix=' + str(prefix), '--disable-tlsopt',
        '--disable-openmp', 'CC=' + a.cc, 'CXX=' + a.cxx, 'FC=' + a.fc, 'F77=' + a.fc,
        'FFLAGS=-O2 -fallow-argument-mismatch', 'LDFLAGS=' + ldflags], source)
    run('dyturbo-build', ['make', '-j' + str(a.jobs)], source)
    run('dyturbo-install', ['make', 'install'], source)
    client_build = a.client_build.resolve()
    run('client-configure', ['cmake', '-S', ROOT / 'adapters/dyturbo', '-B', client_build,
        '-DDYTURBO_ROOT=' + str(prefix), '-DCMAKE_CXX_COMPILER=' + a.cxx], ROOT)
    run('client-build', ['cmake', '--build', client_build, '-j', str(a.jobs)], ROOT)
    report = dict(schema='leptonff-dyturbo-build-1', backend='1.4.2',
        archives={k: dict(url=v[0], sha256=v[1]) for k, v in ARCHIVES.items()},
        prefix=str(prefix), client=str(client_build / 'leptonff_dyturbo'),
        modified_sources={str(x.relative_to(source)): digest(x) for x in (matrix, ct)},
        fast_math=False, openmp=False, commands=commands,
        licenses='DYTurbo GPLv3; external dependencies retain their own licenses; not vendored in leptonff')
    (work / 'build_manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    print('Native client: ' + report['client'])
    print('Install the requested LHAPDF grids separately and set --pdf-data when running.')


if __name__ == '__main__':
    main()
