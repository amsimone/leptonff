#!/usr/bin/env python3
"""Relink an isolated MCFM adapter using an existing compatible GCC build.

Reads host objects/libraries without modifying the installation. Writes all
replacement objects and the executable beneath --output. The core must have
been built with the same libstdc++ ABI (build-host in the reproduction recipe).
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('host_build',type=Path)
    parser.add_argument('--core-build',type=Path,default=Path('build-host'))
    parser.add_argument('--output',type=Path,default=Path('build-mcfm-adapter'))
    parser.add_argument('--photos-build',type=Path,help='Optional local PHOTOS HEPEVT build with a source/ snapshot')
    parser.add_argument('--iteration-audit',action='store_true',help='Record raw rate estimates before histogram iteration weighting')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[2]
    host=args.host_build.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    core=args.core_build.resolve()/'libleptonff.a'
    subprocess.run(['cmake','--build',str(core.parent),'--target','leptonff','--parallel','4'],check=True)
    flags=(host/'CMakeFiles/objlib.dir/flags.make').read_text().splitlines()
    values={line.split(' = ',1)[0]:shlex.split(line.split(' = ',1)[1]) for line in flags if ' = ' in line}
    link=shlex.split((host/'CMakeFiles/mcfm.dir/link.txt').read_text())
    compiler=link[0]
    cxx=compiler.replace('gfortran','g++')
    includes=[f'-I{host}/include',*values['Fortran_INCLUDES']]
    fflags=values['Fortran_FLAGS']
    fflags=[flag for flag in fflags if not flag.startswith('-J')]
    subprocess.run([compiler,*values['Fortran_DEFINES'],*includes,*fflags,f'-J{out}',
        '-c',str(root/'adapters/mcfm/covariance.f90'),'-o',str(out/'covariance.o')],check=True)
    subprocess.run([compiler,*values['Fortran_DEFINES'],*includes,*fflags,f'-J{out}',
        '-c',str(root/'adapters/mcfm/nplotter_user.f90'),'-o',str(out/'nplotter_user.o')],check=True)
    storage_replacement=None
    real_replacement=None
    if args.iteration_audit:
        storage=host.parent/'src/Mods/mod_MCFMStorage.f90'
        source=storage.read_text()
        target='                call accumulateIteration(info%callsPerIt)'
        if source.count(target)!=1:raise RuntimeError('Unexpected host iteration accumulation site')
        source=source.replace(target,'                call leptonff_iteration_audit(stage, currentPart, currentIps, &\n'+
                              '                    info%lastIter+1, info%callsPerIt)\n'+target)
        event_site='                f = wgt * fxn(x,wgt)'
        if source.count(event_site)!=1:raise RuntimeError('Unexpected integration event site')
        source=source.replace(event_site,'                if (bin) call leptonff_cov_begin()\n'+event_site+
                              '\n                if (bin) call leptonff_cov_end()')
        (out/'mod_MCFMStorage.f90').write_text(source)
        subprocess.run([compiler,*values['Fortran_DEFINES'],*includes,*fflags,f'-J{out}',
          '-c',str(out/'mod_MCFMStorage.f90'),'-o',str(out/'mod_MCFMStorage.o')],check=True)
        storage_replacement='CMakeFiles/objlib.dir/src/Mods/mod_MCFMStorage.f90.o'
        if link.count(storage_replacement)!=1:raise RuntimeError('Missing host storage object')
        real_source=host.parent/'src/Procdep/realint.f'
        real_text=real_source.read_text()
        commit='         call threadStorageOp(shtmpcommit)'
        if real_text.count(commit)!=1:raise RuntimeError('Unexpected real/dipole commit site')
        real_text=real_text.replace(commit,'         call leptonff_cov_real_commit()\n'+commit)
        (out/'realint.f').write_text(real_text)
        subprocess.run([compiler,*values['Fortran_DEFINES'],*includes,*fflags,f'-J{out}',
          '-c',str(out/'realint.f'),'-o',str(out/'realint.o')],check=True)
        real_replacement='CMakeFiles/objlib.dir/src/Procdep/realint.f.o'
        if link.count(real_replacement)!=1:raise RuntimeError('Missing host real integration object')

    extra=['-DLEPTONFF_WITH_PHOTOS=1'] if args.photos_build else []
    subprocess.run([cxx,'-std=c++17','-O3','-fPIC',f'-I{root}/include',*extra,
        '-c',str(root/'adapters/mcfm/bridge.cpp'),'-o',str(out/'bridge.o')],check=True)
    replacement='CMakeFiles/objlib.dir/src/User/nplotter_user.f90.o'
    if link.count(replacement)!=1:
        raise RuntimeError('Expected exactly one MCFM user-hook object')
    cmd=[]
    for i,token in enumerate(link):
        if token==replacement:cmd.append(str(out/'nplotter_user.o'))
        elif token==storage_replacement:cmd.append(str(out/'mod_MCFMStorage.o'))
        elif token==real_replacement:cmd.append(str(out/'realint.o'))
        elif i and link[i-1]=='-o':cmd.append(str(out/'mcfm_leptonff'))
        elif token.endswith(('.o','.a','.dylib')) and not Path(token).is_absolute():
            cmd.append(str(host/token))
        else:cmd.append(token)
    cmd.extend([str(out/'bridge.o'),str(out/'covariance.o'),str(core)])
    photos_files=[]
    if args.photos_build:
        photos=args.photos_build.resolve()
        includes=[f'-I{p}' for p in sorted((photos/'source').iterdir()) if p.is_dir()]
        subprocess.run([cxx,'-std=c++17','-O3',f'-I{root}/include',*includes,
          '-c',str(root/'adapters/mcfm/photos_comparison.cpp'),'-o',str(out/'photos_comparison.o')],check=True)
        objects=sorted(p for p in photos.glob('*.o') if p.name!='photos_decay.o')
        if not objects:raise RuntimeError('No PHOTOS library objects found')
        cmd.extend([str(out/'photos_comparison.o'),*[str(p) for p in objects]])
        photos_files=objects+sorted((photos/'source').rglob('*.h'))+sorted((photos/'source').rglob('*.cxx'))
    subprocess.run(cmd,cwd=out,check=True)
    if args.photos_build:
        subprocess.run([cxx,'-std=c++17','-O2',str(root/'adapters/mcfm/test_photos.cpp'),
          str(out/'bridge.o'),str(out/'photos_comparison.o'),str(core),*[str(p) for p in objects],
          '-o',str(out/'test_photos')],check=True)
        for process in [1,6,31]:
            subprocess.run([str(out/'test_photos'),str(process)],env={**__import__('os').environ,'OMP_NUM_THREADS':'1'},check=True)
    manifest={'host_build':str(host),'core_library':str(core),
      'host_executable_sha256':digest(host/'mcfm'),
      'host_link_recipe_sha256':digest(host/'CMakeFiles/mcfm.dir/link.txt'),
      'core_library_sha256':digest(core),
      'bridge_sha256':digest(root/'adapters/mcfm/bridge.cpp'),
      'fortran_hook_sha256':digest(root/'adapters/mcfm/nplotter_user.f90'),
      'adapter_executable_sha256':digest(out/'mcfm_leptonff'),
      'compiler':compiler,'method':'read-only reuse of local host objects; replacement analysis hook',
      'accuracy':'FF-only default; optional photonic LL, unitary-conditioned photonic NLL or coupled MSbar NLL spin-correlated hard matching; one-family/variable-flavour QED with finite thresholds; optional LO QCD gluon/quark mixing and sourced HVP anchor with explicit pointlike boundary; mixed kernels absent'}
    manifest['capabilities']=['photonic-nll-unitary-v1','production-archive-v1','common-decay-api-v1']
    if args.photos_build:
        subprocess.run([cxx,'-std=c++17','-O2',str(root/'adapters/mcfm/test_photonic_nll.cpp'),
          str(out/'bridge.o'),str(out/'photos_comparison.o'),str(core),*[str(p) for p in objects],
          '-o',str(out/'test_photonic_nll')],check=True)
    manifest['iteration_audit']=args.iteration_audit
    if args.iteration_audit:
        manifest['iteration_audit_original_sha256']=digest(storage)
        manifest['iteration_audit_source_sha256']=digest(out/'mod_MCFMStorage.f90')
        manifest['real_covariance_original_sha256']=digest(real_source)
        manifest['real_covariance_source_sha256']=digest(out/'realint.f')
    manifest['event_covariance_sha256']=digest(root/'adapters/mcfm/covariance.f90')
    manifest['photos_sha256']={str(p):digest(p) for p in photos_files}
    sys.path.insert(0,str(root/'python'))
    from check_reference import source_digest
    manifest['source_sha256']=source_digest(root)
    manifest['photos_comparison_sha256']=digest(root/'adapters/mcfm/photos_comparison.cpp') if args.photos_build else None
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(out/'mcfm_leptonff')

if __name__=='__main__':main()
