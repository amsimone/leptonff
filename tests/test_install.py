"""Build C and Fortran consumers against a relocated installed package."""
import argparse
import os
from pathlib import Path
import shutil
import shlex
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('build', type=Path)
p.add_argument('--cxx', required=True)
p.add_argument('--fortran')
p.add_argument('--tools', default='ON')
p.add_argument('--core-only', action='store_true')
a = p.parse_args()
root = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='leptonff-install-') as temp:
    base = Path(temp)
    def run(command):
        result = subprocess.run(list(map(str, command)), capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(result.stdout+'\n'+result.stderr)
        return result.stdout
    run(['cmake', '--install', a.build.resolve(), '--prefix', base/'installed'])
    shutil.move(base/'installed', base/'relocated')
    examples = ['fragmentation','convolution']
    if not a.core_only:
        examples += ['c_api','cpp','response_table']+(['fortran'] if a.fortran else [])
    for example in examples:
        command = ['cmake', '-S', base/'relocated/share/doc/leptonff/examples'/example, '-B', base/example,
                   '-DCMAKE_PREFIX_PATH='+str(base/'relocated'), '-DCMAKE_CXX_COMPILER='+a.cxx]
        if example == 'fortran':
            command.append('-DCMAKE_Fortran_COMPILER='+a.fortran)
        run(command)
        run(['cmake', '--build', base/example, '--parallel', '2'])
        executable = {'fragmentation':'leptonff_fragmentation_example',
                      'convolution':'leptonff_convolution_example',
                      'c_api':'decay_example','cpp':'decay_cpp','fortran':'decay_fortran','response_table':'response_example'}[example]
        output = run([base/example/executable])
        expected = {'fragmentation':'fragmentation moment','convolution':'convolution','response_table':'matched contraction'}.get(example,'Born acceptance')
        if expected not in output:
            raise AssertionError('Installed consumer did not evaluate the documented example')
        print(example, output.strip())
    if not (base/'relocated/share/doc/leptonff/benchmarks/canonical/order1.json').is_file():
        raise AssertionError('Canonical reference data not installed')
    if shutil.which('pkg-config'):
        pc=next((base/'relocated').glob('lib*/pkgconfig/leptonff.pc'))
        environment={**os.environ,'PKG_CONFIG_PATH':str(pc.parent)}
        flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','leptonff'],
            env=environment,text=True))
        client=base/'pkgconfig.cpp'
        client.write_text('#include <leptonff/core.hpp>\n#include <iostream>\nint main(){std::cout<<leptonff::version();}\n')
        run([a.cxx,client,'-std=c++17',*flags,'-Wl,-rpath,'+str(pc.parent.parent),'-o',base/'pkgconfig'])
        expected_version=subprocess.check_output(['pkg-config','--modversion','leptonff'],env=environment,text=True).strip()
        if run([base/'pkgconfig']).strip()!=expected_version:
            raise AssertionError('pkg-config linked the wrong native version')
    if a.core_only:
        probe=base/'missing-component'
        probe.mkdir()
        (probe/'CMakeLists.txt').write_text('cmake_minimum_required(VERSION 3.16)\nproject(probe LANGUAGES NONE)\nfind_package(leptonff 0.27 CONFIG REQUIRED COMPONENTS matched)\n')
        result=subprocess.run(['cmake','-S',str(probe),'-B',str(probe/'build'),
            '-DCMAKE_PREFIX_PATH='+str(base/'relocated')],capture_output=True,text=True)
        if result.returncode==0:
            raise AssertionError('Core-only installation claimed the matched component')
    if a.tools == 'ON':
        if not (base/'relocated/bin/leptonff_replay').is_file():
            raise AssertionError('Replay executable not installed')
        if 'Usage:' not in run([base/'relocated/bin/leptonff_replay','--help']):
            raise AssertionError('Installed replay runtime is not relocatable')
print('Relocated package works independently of source/build include paths')
