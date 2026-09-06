
import sys, time
tA = time.perf_counter()
from Bio.SeqUtils import MeltingTemp as mt
tB = time.perf_counter()
seqs = [l.strip() for l in open(sys.argv[1]) if l.strip()]
t0 = time.perf_counter()
tm = [mt.Tm_NN(s, nn_table=mt.DNA_NN4, saltcorr=3,
               Na=50, K=0, Tris=0, Mg=0, dNTPs=0,
               dnac1=25, dnac2=25, selfcomp=False) for s in seqs]
t1 = time.perf_counter()
open(sys.argv[2], "w").write("\n".join("%.10f" % v for v in tm) + "\n")
t2 = time.perf_counter()
tim = ("LOAD_SECONDS %.4f\nREAD_SECONDS %.4f\n"
       "COMPUTE_SECONDS %.4f\nWRITE_SECONDS %.4f\n"
       % (tB - tA, t0 - tB, t1 - t0, t2 - t1))
if len(sys.argv) > 3:
    open(sys.argv[3], "w").write(tim)
print(tim)

