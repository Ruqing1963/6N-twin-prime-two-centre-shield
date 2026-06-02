# Arithmetic lattice interferometer (S10, memory-efficient).
# Directly MEASURE right-centre survival P(N+d twin | N twin, group), grouping the
# left twin centre N by divisibility by {5,7}, restricted to a high-omega band,
# and compare to the CRT independence baseline (= same probability over the whole
# band). No "independent product" model assumed. Reveals the real two-centre
# coupling. (S9 finding to confirm: the coupling on 42 is 5-driven, not 7-driven.)
#
# Memory-safe: stores twin centres as a sorted numpy array + parallel arrays for
# omega and divisibility, and tests "N+d is twin" by binary search. ~ a few hundred
# MB for S10 rather than multi-GB dicts.
import numpy as np, math, os, time
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0].astype(np.int64)
MAXK=int(os.environ.get("MAXK",10))
LO=10**(MAXK-1)//6+1; HI=10**MAXK//6; SEG=4_000_000
PB=int(math.isqrt(6*HI+250))+1; BP=primes_upto(PB)
GAPS=[7,35]  # 42, 210

N_list=[]; om_list=[]
n=LO; t0=time.time()
while n<=HI:
    nh=min(n+SEG,HI+1); sz=nh-n
    rem=np.arange(n,nh,dtype=np.int64); ob=np.zeros(sz,np.int16)
    for p in BP:
        if p*p>nh-1: break
        f=((n+p-1)//p)*p
        if f>=nh: continue
        idx=np.arange(f-n,sz,p)
        if idx.size==0: continue
        sub=rem[idx]; m=(sub%p)==0
        while m.any(): sub[m]//=p; m=(sub%p)==0
        rem[idx]=sub
        if p>3: ob[idx]+=1
    ob[rem>1]+=1
    Narr=np.arange(n,nh,dtype=np.int64)
    vlo=6*n-1; vhi=6*(nh-1)+1; span=vhi-vlo+1
    comp=np.zeros(span,bool); sq=int(math.isqrt(vhi))+1
    for p in BP:
        if p>sq: break
        st=max(p*p,((vlo+p-1)//p)*p)
        if st>vhi: continue
        comp[st-vlo:span:p]=True
    tw=(~comp[(6*Narr-1)-vlo])&(~comp[(6*Narr+1)-vlo])
    pos=np.nonzero(tw)[0]
    N_list.append(Narr[pos]); om_list.append(ob[pos])
    n=nh
N_arr=np.concatenate(N_list)          # sorted ascending by construction
om_arr=np.concatenate(om_list).astype(np.int16)
print(f"S{MAXK} twins {len(N_arr):,}; scan {time.time()-t0:.0f}s")

def is_twin(vals):
    # vectorised membership test against the sorted N_arr
    idx=np.searchsorted(N_arr, vals)
    idx=np.clip(idx,0,len(N_arr)-1)
    return N_arr[idx]==vals

div5=(N_arr%5==0); div7=(N_arr%7==0)
band=(om_arr>=4)&(om_arr<=6)
print(f"high-omega band (omega 4..6): {band.sum():,} twin centres")

groups={
 'A(5,!7)': div5&~div7&band,
 'B(7,!5)': div7&~div5&band,
 'C(5&7)' : div5&div7&band,
 'D(!5,!7)':~div5&~div7&band,
}
for d in GAPS:
    sixd=6*d
    Nb=N_arr[band]
    right_twin=is_twin(Nb+d)
    base=right_twin.mean()
    print(f"\n=== 6dN={sixd} (d={d}) ; baseline P(N+d twin | N twin, band) = {base:.4f} ===")
    print(f"  {'group':>10}{'n_left':>11}{'P(right twin)':>15}{'vs baseline':>13}")
    for g,maskg in groups.items():
        Ng=N_arr[maskg]
        if len(Ng)==0: continue
        p=is_twin(Ng+d).mean()
        print(f"  {g:>10}{len(Ng):>11,}{p:>15.4f}{p/base:>13.3f}")
