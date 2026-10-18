#!/usr/bin/env python3
"""Build a ROOT-free FSR reference from checksum-pinned ReneSANCe modules.

External/generated GPL sources stay in the ignored build/source directories.
Only the published Born, virtual/soft, and hard matrix elements are used;
the benchmark supplies integration, histograms, and independent QED decays.
"""
import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.request

HERE = Path(__file__).resolve().parent
SHA256 = '260358cbc558fa3edac0d4e2f5447a5ab6296396100abaf38b6b83a7597b7165'
URL = 'https://renesance.hepforge.org/downloads?f=ReneSANCe-1.3.11.tar.gz'


def defaults(source):
    values = {}
    for name in ('renesance.schema', 'model.schema', 'foam.schema'):
        for key, body in re.findall(r'^\s*(\w+)\s*:\s*\{([^{}]*)\}',
                                    (source/'Generator/schema'/name).read_text(), re.M):
            match = re.search(r'\bdefault\s*:\s*([^,}]+)', body)
            if match:
                values[key] = match[1].strip().strip('"')
    return values


def neutral_fsr_integrated(text, channel):
    """Select the final-state terms in the release's analytic dσ/ds'.

    Unlike the differential nc_ha routine, nc_ha_*_1spr ignores iqed.
    v0s/vas multiply the current at s (FSR); v0spr/vaspr describe ISR,
    and a1sspr/a2sspr its interference. Preserve every FSR coefficient.
    """
    name = 'nc_ha_'+channel+'_1spr'
    block = text[text.index('      subroutine '+name):]
    start = block.index('        hard=\n')
    end = block.index('\n        hard = hard*conhc*cfprime', start)
    expression = block[start+len('        hard=\n'):end]
    terms = re.split(r'(?m)(?=^     &   \+(?:v0s|vas|a1sspr|a2sspr)\b|^     &   \+(?:v0spr|vaspr)\b)', expression)
    retained = [t for t in terms if re.match(r'     &   \+(?:v0s|vas)\*', t)]
    if len(retained)!=8 or len([t for t in terms if t.strip()])!=14:
        raise ValueError('Pinned NC integrated-expression structure changed')
    block = block[:start]+'        hard=\n'+''.join(retained)+block[end:]
    notice='* Derived from ReneSANCe 1.3.11 / SANC (GPLv3).\n* Copyright SANC Project Team, 2002-2012.\n* FSR terms selected without changing their coefficients.\n'
    return notice+block.replace('subroutine '+name, 'subroutine nc_fsr_ha_'+channel+'_1spr', 1)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive', type=Path)
    p.add_argument('--build', type=Path, default=HERE.parents[1]/'build-renesance')
    p.add_argument('--lhapdf', type=Path, required=True)
    p.add_argument('--looptools', type=Path,
                   help='Reuse a built LoopTools directory; otherwise build the bundled release')
    p.add_argument('--leptonff', type=Path, required=True, help='Installed leptonff prefix')
    p.add_argument('--photos-source', type=Path, required=True)
    p.add_argument('--photos-library', type=Path, required=True)
    p.add_argument('--jobs', type=int, default=4)
    args = p.parse_args()
    build = args.build.resolve(); build.mkdir(parents=True, exist_ok=True)
    archive = args.archive or build/'ReneSANCe-1.3.11.tar.gz'
    if not archive.exists():
        with urllib.request.urlopen(URL, timeout=60) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != SHA256:
            raise ValueError('ReneSANCe download checksum mismatch')
        archive.write_bytes(data)
    if hashlib.sha256(archive.read_bytes()).hexdigest() != SHA256:
        raise ValueError('ReneSANCe archive checksum mismatch')
    source = build/'ReneSANCe-1.3.11'
    if not source.exists():
        with tarfile.open(archive) as tar:
            for member in tar.getmembers():
                path = Path(member.name)
                if path.is_absolute() or '..' in path.parts or member.issym() or member.islnk():
                    raise ValueError('Unsafe source archive member')
            tar.extractall(build)
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            if member.isfile() and '/Generator/' in member.name:
                original=tar.extractfile(member).read()
                if (build/member.name).read_bytes()!=original:
                    raise ValueError('Extracted ReneSANCe input changed: '+member.name)
    generator = source/'Generator'
    # These two ROOT conveniences are not used in the FSR-only path. Fail if
    # an unsupported ISR subtraction calls them, rather than silently faking it.
    (build/'TMath.h').write_text('#pragma once\n#include <stdexcept>\n'
        'namespace TMath { inline double DiLog(double) {throw std::runtime_error('
        '"ISR subtraction is unsupported by the FSR-only driver");} }\n')
    # Select the hadron Process declaration in the upstream initialisation unit.
    init = (generator/'src/ModulesInterface.cpp').read_text()
    (build/'ModulesInterface.cpp').write_text(init.replace('#include "Process.h"',
                                                          '#include "Process_pp.h"'))
    values = defaults(source)
    (build/'defaults.hpp').write_text('#pragma once\n#include <map>\n#include <string>\n'
        'inline std::map<std::string,std::string> sanc_defaults() { return {\n'+
        ',\n'.join('{'+json.dumps(k)+','+json.dumps(v)+'}' for k,v in sorted(values.items()))+
        '\n}; }\n')
    fc = shutil.which('gfortran'); cxx = shutil.which('clang++') or shutil.which('c++')
    if not fc or not cxx:
        raise RuntimeError('gfortran and a C++ compiler are required')
    if args.looptools:
        lt=args.looptools.resolve()
    else:
        lt=build/'LoopTools'
        if not lt.exists():
            with tarfile.open(source/'external/LoopTools.tar.gz') as tar:
                for member in tar.getmembers():
                    path=Path(member.name)
                    if path.is_absolute() or '..' in path.parts or member.issym() or member.islnk():
                        raise ValueError('Unsafe LoopTools archive member')
                tar.extractall(build)
        if not (lt/'libooptools.a').exists():
            with (build/'looptools-build.log').open('w') as log:
                subprocess.run(['./configure'],cwd=lt,stdout=log,stderr=subprocess.STDOUT,check=True)
                subprocess.run(['make',f'-j{args.jobs}','libooptools.a',
                    'FFLAGS=-O2 -fPIC -std=legacy -fallow-argument-mismatch -DQEXT=DBLE',
                    'CFLAGS=-O2 -fPIC -DHAVE_UNDERSCORE'],cwd=lt,stdout=log,stderr=subprocess.STDOUT,check=True)
    pdf=args.lhapdf.resolve()
    ff = [fc, '-O2', '-fPIC', '-std=legacy', '-fallow-argument-mismatch',
          '-ffixed-line-length-none', '-fno-range-check', '-I'+str(generator/'SANC'),
          '-I'+str(lt/'include')]
    fortran = sorted((generator/'SANC').glob('*.f'))+sorted((generator/'SANC').glob('*.F'))
    fortran += sorted((generator/'modules_pp').glob('*.f'))+sorted((generator/'modules_pp').glob('*.F'))
    for channel in ('1313_1616','1414_1616'):
        original=generator/'modules_pp'/('nc_br_'+channel+'.F')
        generated=build/('nc_fsr_'+channel+'.F')
        generated.write_text(neutral_fsr_integrated(original.read_text(),channel))
        fortran.append(generated)
    commands = []
    objects = []
    for src in fortran:
        obj = build/(src.parent.name+'_'+src.name+'.o'); objects.append(obj)
        if not obj.exists() or obj.stat().st_mtime < src.stat().st_mtime:
            commands.append(ff+['-c',str(src),'-o',str(obj)])
    def compile_one(command):
        r = subprocess.run(command, capture_output=True, text=True)
        return command, r.returncode, r.stdout+r.stderr
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(compile_one, commands))
    (build/'fortran-build.log').write_text('\n'.join(' '.join(c)+'\n'+s for c,_,s in results))
    if any(code for _,code,_ in results):
        raise RuntimeError('Fortran compilation failed; see fortran-build.log')
    lib = build/'libsanc_reference.a'
    lib.unlink(missing_ok=True)
    subprocess.run(['ar','rcs',str(lib),*[str(o) for o in objects]],check=True)
    cpp = [generator/'src'/n for n in ('Info.cpp','Process_pp.cpp','Vec4.cpp')]
    cpp += [build/'ModulesInterface.cpp', HERE/'reference.cpp', HERE/'pythia_decay.cpp',
            HERE.parent/'dyturbo/photos_decay.cpp']
    cflags = [cxx,'-O2','-std=c++17','-fPIC','-I'+str(build),'-iquote'+str(generator/'include'),
              '-I'+str(pdf/'include'),'-I'+str(args.leptonff.resolve()/'include'),
              '-I'+str(HERE.parent/'dyturbo')]
    photos = args.photos_source.resolve()
    cflags += ['-I'+str(d) for d in [photos,*photos.iterdir()] if d.is_dir()]
    cflags += [f for f in subprocess.check_output(['pythia8-config','--cxxflags'],text=True).split()
               if not f.startswith('-std=')]
    cpp_commands = []
    for src in cpp:
        obj = build/(src.stem+'.o'); objects.append(obj)
        cpp_commands.append(cflags+['-c',str(src),'-o',str(obj)])
    results = [compile_one(c) for c in cpp_commands]
    (build/'cpp-build.log').write_text('\n'.join(' '.join(c)+'\n'+s for c,_,s in results))
    if any(code for _,code,_ in results):
        raise RuntimeError('C++ compilation failed; see cpp-build.log')
    runtime = Path(subprocess.check_output([fc,'-print-file-name='+('libgfortran.dylib' if sys.platform=='darwin' else 'libgfortran.so')],text=True).strip()).parent
    command = [cxx,'-o',str(build/'leptonff_fsr_reference'),
               *[str(build/(src.stem+'.o')) for src in cpp], str(lib),
               str(lt/'libooptools.a'),
               str(args.photos_library.resolve()),'-L'+str(pdf/'lib'),'-lLHAPDF',
               '-Wl,-rpath,'+str(pdf/'lib'),'-L'+str(runtime),'-lgfortran',
               '-Wl,-rpath,'+str(runtime),'-L'+str(args.leptonff.resolve()/'lib'),
               '-lleptonff','-Wl,-rpath,'+str(args.leptonff.resolve()/'lib')]
    command += subprocess.check_output(['pythia8-config','--libs'],text=True).split()
    r = subprocess.run(command,capture_output=True,text=True)
    (build/'link.log').write_text(r.stdout+r.stderr)
    if r.returncode:
        raise RuntimeError('Link failed; see link.log')
    (build/'manifest.json').write_text(json.dumps({'archive_sha256':SHA256,'url':URL,
        'version':'1.3.11','source_sha256':{str(s.relative_to(source)):
        hashlib.sha256(s.read_bytes()).hexdigest() for s in fortran+cpp[:3] if s.is_relative_to(source)},
        'generated_sources':{s.name:hashlib.sha256(s.read_bytes()).hexdigest() for s in fortran if not s.is_relative_to(source)},
        'link_command':command,'looptools_library_sha256':hashlib.sha256((lt/'libooptools.a').read_bytes()).hexdigest(),
        'photos_library_sha256':hashlib.sha256(args.photos_library.read_bytes()).hexdigest(),
        'adapter_sources':{s.name:hashlib.sha256(s.read_bytes()).hexdigest() for s in [HERE/'build.py',*cpp[4:],
            HERE/'pythia_decay.hpp',HERE.parent/'dyturbo/mass_radiator.hpp',HERE.parent/'dyturbo/photos_decay.hpp']}},indent=2)+'\n')
    print(build/'leptonff_fsr_reference')


if __name__ == '__main__':
    main()
