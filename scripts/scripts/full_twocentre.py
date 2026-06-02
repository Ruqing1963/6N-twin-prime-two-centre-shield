#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Two-centre shield law: modular-shift gain x cross-prime hedge  (Part V)
================================================================================
For a centre-step d and prime q>3, the two-centre shield is
    S(d,q) = P(N+d twin | q|N, N twin) / P(N+d twin | N twin).
This script measures S(d,q) over a high-omega band of S_K and decomposes it as
    S(d,q) = g(d,q) * h(d,q),  where
  g (modular-shift gain): q|N forces N+d ≡ d (mod q); if d is q-safe the right
      centre is certainly q-safe, vs its MEASURED baseline q-safe rate beta_q(d).
      g = 1{d not in dead(q)} / beta_q(d).
  h (cross-prime hedge): conditioning on q|N moves the right centre's residues at
      the OTHER primes q'; h = prod_{q'!=q} [ P(N+d q'-safe | q|N) / beta_{q'}(d) ].
Both factors use the MEASURED modular distribution of twin centres (which carries
the Part-I enrichment); an idealised uniform baseline overstates g and misses h.

Result: g*h reproduces the measured shield for 6dN = 42 and 210 to <=1.5% in S10
(<=3% in S9). The hedge is <1 for 42 (q=5 collision) and >1 for 210. The 42
shield is a q=5 modular-shift effect, NOT q=7 (the gap's own factor is inert).

SCOPE: this is the right-centre conditional survival law. Bridging it to the
omega-stratified gap preference r(d|omega) of Part IV (and the fate of the 42
residual) is left open in the paper.

USAGE
  python full_twocentre.py            # default S10 (~15-20 min)
  MAXK=9 python full_twocentre.py     # S9 (faster, for validation)
Requires: numpy. Emits tc_plot_S{K}.csv.
================================================================================
"""
# Full two-centre shield model = (own modular-shift gain) x (hedge: product over
# other primes of M's safety-rate ratio). Predict the full shield for every q|N
# and compare to measurement, for d=7 and d=35. If it closes, we have the
# quantitative two-centre law.
import numpy as np, math, os
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0].astype(np.int64)
MAXK=int(os.environ.get("MAXK",10))
LO=10**(MAXK-1)//6+1; HI=10**MAXK//6; SEG=4_000_000
PB=int(math.isqrt(6*HI+250))+1; BP=primes_upto(PB)
N_list=[]; om_list=[]
n=LO
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
N_arr=np.concatenate(N_list); om_arr=np.concatenate(om_list).astype(np.int16)
print(f"S{MAXK} twins {len(N_arr):,}")
def is_twin(vals):
    idx=np.searchsorted(N_arr,vals); idx=np.clip(idx,0,len(N_arr)-1)
    return N_arr[idx]==vals
band=(om_arr>=4)&(om_arr<=6); Nb=N_arr[band]
def dead_set(q): inv=pow(6,-1,q); return list({inv%q,(-inv)%q})
QS=[5,7,11,13,17,19,23,29,31,37,41,43]
for d in [7,35]:
    M_all=Nb+d
    base_safe={q:(~np.isin(M_all%q,dead_set(q))).mean() for q in QS}
    base=is_twin(Nb+d).mean()
    print(f"\n=== d={d} (gap {6*d}) full two-centre model vs measurement ===")
    print(f"  {'q':>4}{'own gain':>10}{'hedge':>9}{'predicted':>11}{'measured':>10}{'err%':>7}")
    maxerr=0
    for q in QS:
        selq=Nb[Nb%q==0]
        if len(selq)<1500: continue
        Mq=selq+d
        own=(~np.isin(Mq%q,dead_set(q))).mean()/base_safe[q] if base_safe[q]>0 else 0
        hedge=1.0
        for qp in QS:
            if qp==q: continue
            sc=(~np.isin(Mq%qp,dead_set(qp))).mean()
            hedge*= (sc/base_safe[qp]) if base_safe[qp]>0 else 1
        pred=own*hedge
        meas=is_twin(selq+d).mean()/base
        err=100*(pred-meas)/meas if meas>0 else 0
        if abs(err)>abs(maxerr): maxerr=err
        print(f"  {q:>4}{own:>10.3f}{hedge:>9.3f}{pred:>11.3f}{meas:>10.3f}{err:>7.1f}")
    print(f"  max |error| = {abs(maxerr):.1f}%")


# ---- emit CSV for plotting ----
import csv as _csv
with open(f'tc_plot_S{MAXK}.csv','w',newline='') as _f:
    _w=_csv.writer(_f); _w.writerow(['gap','q','own_gain','hedge','predicted','measured'])
    for d in [7,35]:
        M_all=Nb+d
        base_safe={q:(~np.isin(M_all%q,dead_set(q))).mean() for q in QS}
        base=is_twin(Nb+d).mean()
        for q in QS:
            selq=Nb[Nb%q==0]
            if len(selq)<1500: continue
            Mq=selq+d
            own=(~np.isin(Mq%q,dead_set(q))).mean()/base_safe[q] if base_safe[q]>0 else 0
            hedge=1.0
            for qp in QS:
                if qp==q: continue
                sc=(~np.isin(Mq%qp,dead_set(qp))).mean()
                hedge*=(sc/base_safe[qp]) if base_safe[qp]>0 else 1
            pred=own*hedge; meas=is_twin(selq+d).mean()/base
            if pred>0:
                _w.writerow([6*d,q,f'{own:.4f}',f'{hedge:.4f}',f'{pred:.4f}',f'{meas:.4f}'])
print(f"\n[ok] wrote tc_plot_S{MAXK}.csv")
