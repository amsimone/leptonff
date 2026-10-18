"""Run after pip install . and setting LEPTONFF_LIBRARY to a shared library."""
from leptonff import Fragmentation

with Fragmentation() as ff:
    print("N=2 fragmentation moment", ff.moment(2))
    print("Endpoint bin [0.9,1]", ff.integral(.9, 1))
    print("Convolution at x=0.5", ff.convolve(lambda x: x * (1 - x), .5))
