#!/usr/bin/env python3
"""Read-only comparison with an installed Jegerlehner hadr5x19 routine.

Compiles the supplied external file without modifying or copying it; records
its hash and exports a separate alternative anchor. The table interpolation
in this routine is single precision, including its double-argument wrapper.
"""
import argparse,hashlib,json,math,subprocess
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source',type=Path)
    p.add_argument('--compiler',default='/opt/homebrew/bin/gfortran-16')
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    driver=out/'probe.f90';driver.write_text('''program probe
implicit none
real(8) :: e,st,der,stat,syst,deg,gs,gy
integer :: i
real(8),parameter :: scales(4)=[91.1876d0,80.385d0,-2d0,-4d0]
st=0.23153d0
do i=1,4
e=scales(i)
call dhadr5x19(e,st,der,stat,syst,deg,gs,gy)
write(*,'(4es25.16)') e,der,stat,syst
end do
end program
''')
    exe=out/'probe';cmd=[a.compiler,'-O2','-ffixed-line-length-none','-fallow-argument-mismatch',str(a.source.resolve()),str(driver),'-o',str(exe)]
    subprocess.run(cmd,check=True,capture_output=True,text=True)
    raw=subprocess.run([str(exe)],check=True,capture_output=True,text=True).stdout
    rows=[dict(zip(['energy_GeV','delta_had5_os','statistical_error','systematic_error'],map(float,line.split()))) for line in raw.splitlines()]
    anchor=rows[0];error=math.hypot(anchor['statistical_error'],anchor['systematic_error'])
    # Different sourced evaluations need agreement within their uncertainties,
    # not bitwise equality. PDG 2025 fit input: .02760 +/- .00010.
    pull=abs(anchor['delta_had5_os']-.02760)/math.hypot(error,.00010)
    if len(rows)!=4 or pull>3 or any(not 0<r['delta_had5_os']<.1 for r in rows):raise RuntimeError('External HVP reference failure')
    report={'passed':True,'source':str(a.source.resolve()),'source_sha256':hashlib.sha256(a.source.read_bytes()).hexdigest(),
        'compile_command':cmd,'points':rows,'pdg2025_unpaired_input_pull':pull,
        'scope':'Alternative Jegerlehner 2019 OS HVP evaluation; no fitted real-fragmentation boundary or complete EW coupling'}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    (out/'jegerlehner19.dat').write_text(f'# Generated from installed hadr5x19; source sha256 {report["source_sha256"]}\nreference_scale 91.1876\ndelta_had5_os {anchor["delta_had5_os"]:.17g}\nuncertainty {error:.17g}\nalpha_s_reference 0.1193\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
