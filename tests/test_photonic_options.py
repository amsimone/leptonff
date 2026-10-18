"""Reject unsupported or inconsistent photonic-NLL launches before any run."""
from pathlib import Path
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[1]
runner = root/'adapters/mcfm/run_reference.py'
with tempfile.TemporaryDirectory() as directory:
    base = Path(directory)
    cases = [(['--photonic-nll'], 'requires the adapter'),
             (['--photonic-nll', '--photonic-ll'], 'not allowed with argument'),
             (['--photonic-nll', '--dress-hard', '--hard-samples', '4', '--native'], 'requires the adapter'),
             (['--photonic-nll', '--dress-hard', '--hard-samples', '4', '--coupled-hard'], 'requires coupled evolution'),
             (['--photonic-nll', '--dress-hard', '--hard-samples', '4', '--alpha-scheme', 'on-shell'], 'fixed alpha'),
             (['--photonic-nll', '--dress-hard', '--hard-samples', '4'], 'photonic-nll-unitary-v1')]
    for index, (options, message) in enumerate(cases):
        output = base/f'run{index}'
        command = [sys.executable, str(runner), str(base/'no-host'), 'wplus',
                   '--output', str(output), '--adapter', str(base/'old-adapter'), *options]
        run = subprocess.run(command, text=True, capture_output=True)
        assert run.returncode != 0 and message in run.stderr, run.stderr
        assert not output.exists(), 'Invalid theory mode created a run directory'
print('Six incompatible/unsupported NLL launches rejected before creating output or starting the host')
