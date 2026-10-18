"""Check distributed canonical numbers against the numerical library."""
import json
from pathlib import Path
import subprocess
import sys
reference=json.loads((Path(__file__).resolve().parents[1]/'benchmarks/canonical/order1.json').read_text())
actual=json.loads(subprocess.check_output([sys.argv[1]],text=True))
assert actual.keys()==reference['values'].keys()
for key,value in reference['values'].items():
    assert abs(actual[key]-value)<reference['absolute_tolerance'], (key,actual[key],value)
print('Canonical analytic fragmentation, endpoint and convolution benchmarks passed')
