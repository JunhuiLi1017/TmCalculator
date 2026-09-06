
import sys, time
from Bio.SeqUtils import MeltingTemp as mt
seqs = [l.strip() for l in open(sys.argv[1]) if l.strip()]
t0 = time.perf_counter()
tm = [mt.Tm_NN(s, nn_table=mt.DNA_NN4, saltcorr=0,
               Na=50, K=0, Tris=0, Mg=0, dNTPs=0,
               dnac1=250, dnac2=250, selfcomp=False) for s in seqs]
t1 = time.perf_counter()
open(sys.argv[2], "w").write("\n".join("%.10f" % v for v in tm) + "\n")
print("COMPUTE_SECONDS %.4f" % (t1 - t0))

