#!/usr/bin/env python3
"""Build the official KKMCee 4.32.01 Fortran core without ROOT or EW tables."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
COMMIT='527ed9d6f8c0080a256e77ae0d7adbf4c5207eea'
URL='https://github.com/KrakowHEPSoft/KKMCee.git'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=HERE.parents[1]/'external/kkmcee-f77-audit')
    p.add_argument('--build',type=Path,default=HERE.parents[1]/'build-kkmcee-fsr')
    p.add_argument('--jobs',type=int,default=4)
    a=p.parse_args();source=a.source.resolve();build=a.build.resolve()
    if not source.exists():
        subprocess.run(['git','clone','--depth','1','--branch','v4.32.01',URL,str(source)],check=True)
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip()!=COMMIT:
        raise ValueError('KKMC source commit differs from the pinned release')
    if subprocess.check_output(['git','diff','--name-only','HEAD'],cwd=source,text=True).strip():
        raise ValueError('Tracked KKMC sources have changed')
    build.mkdir(parents=True,exist_ok=True)
    fc=shutil.which('gfortran');cxx=shutil.which('clang++') or shutil.which('c++')
    if not fc or not cxx:raise RuntimeError('Need gfortran and a C++ compiler')
    flags='-O2 -fPIC -std=legacy -fallow-argument-mismatch -fno-second-underscore -fno-automatic -ffixed-line-length-132 -fbackslash'
    with (build/'upstream-build.log').open('w') as log:
        for folder,makefile in [('KK2f','KKMakefile'),('glibk','KKMakefile'),('jetset','KKMakefile'),
                                ('bornv','KKMakefile'),('tauola','makefile'),('photos','KKMakefile')]:
            subprocess.run(['make','-f',makefile,f'-j{a.jobs}','FFLAGS='+flags,'CC='+cxx],
                           cwd=source/folder,stdout=log,stderr=subprocess.STDOUT,check=True)
    runtime=Path(subprocess.check_output([fc,'-print-file-name='+('libgfortran.dylib' if sys.platform=='darwin' else 'libgfortran.so')],text=True).strip()).parent
    libraries=[source/f for f in ['KK2f/KK2f.a','glibk/glibk.a','jetset/jetset.a',
                                 'bornv/bornv-tabs.a','tauola/glib.a','photos/glib.a']]
    # Repeated archives resolve cyclic Fortran utility references on ELF too.
    command=[cxx,'-O2','-std=c++17',str(HERE/'reference.cpp'),'-o',str(build/'kkmcee_fsr_reference'),
             *map(str,libraries if sys.platform=='darwin' else libraries+libraries),'-L'+str(runtime),'-lgfortran',
             '-Wl,-rpath,'+str(runtime)]
    with (build/'link.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    tracked=subprocess.check_output(['git','ls-files'],cwd=source,text=True).splitlines()
    manifest=dict(schema='leptonff-kkmcee-build-v1',url=URL,tag='v4.32.01',commit=COMMIT,
                  source_tree_sha256=hashlib.sha256(b''.join(
                      f.encode()+b'\0'+(source/f).read_bytes() for f in tracked if (source/f).is_file())).hexdigest(),
                  driver_sha256=hashlib.sha256((HERE/'reference.cpp').read_bytes()).hexdigest(),
                  defaults_sha256=hashlib.sha256((source/'.KK2f_defaults').read_bytes()).hexdigest(),
                  command=command,fortran_flags=flags)
    (build/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(build/'kkmcee_fsr_reference')


if __name__=='__main__':main()
