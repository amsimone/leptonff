"""Build C and Fortran consumers against a relocated installed package."""
import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('build', type=Path)
p.add_argument('--cxx', required=True)
p.add_argument('--fortran')
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
    examples = ['c_api','cpp']+(['fortran'] if a.fortran else [])
    for example in examples:
        command = ['cmake', '-S', root/'examples'/example, '-B', base/example,
                   '-DCMAKE_PREFIX_PATH='+str(base/'relocated'), '-DCMAKE_CXX_COMPILER='+a.cxx]
        if example == 'fortran':
            command.append('-DCMAKE_Fortran_COMPILER='+a.fortran)
        run(command)
        run(['cmake', '--build', base/example, '--parallel', '2'])
        executable = {'c_api':'decay_example','cpp':'decay_cpp','fortran':'decay_fortran'}[example]
        output = run([base/example/executable])
        if 'Born acceptance' not in output:
            raise AssertionError('Installed consumer did not evaluate a decay')
        print(example, output.strip())
    if not (base/'relocated/bin/leptonff_replay').is_file():
        raise AssertionError('Replay executable not installed')
    if 'Usage:' not in run([base/'relocated/bin/leptonff_replay','--help']):
        raise AssertionError('Installed replay runtime is not relocatable')
print('Relocated package works independently of source/build include paths')
