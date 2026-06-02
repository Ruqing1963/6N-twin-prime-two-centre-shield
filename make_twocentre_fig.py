#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the two-panel two-centre figure (6dN=42, 210) from ../data/tc_plot_S10.csv
(produced by full_twocentre.py with default MAXK=10). Shows measured shield,
the law g*h, and the two factors (gain, hedge) separately.
"""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows=list(csv.DictReader(open('../data/tc_plot_S10.csv')))
def series(gap):
    q=[];pred=[];meas=[];own=[];hed=[]
    for r in rows:
        if int(r['gap'])==gap:
            q.append(int(r['q'])); pred.append(float(r['predicted'])); meas.append(float(r['measured']))
            own.append(float(r['own_gain'])); hed.append(float(r['hedge']))
    return np.array(q),np.array(pred),np.array(meas),np.array(own),np.array(hed)
fig,axes=plt.subplots(1,2,figsize=(14,5.4))
for ax,gap,col,errlbl in [(axes[0],42,'#185FA5','1.0\\%'),(axes[1],210,'#c0392b','1.5\\%')]:
    q,pred,meas,own,hed=series(gap)
    x=np.arange(len(q))
    ax.plot(x,meas,'o-',color=col,lw=2.2,ms=8,label='measured shield',zorder=4)
    ax.plot(x,pred,'s--',color='#e8845b',lw=1.6,ms=7,label='two-centre law (own$\\times$hedge)',zorder=3)
    ax.plot(x,own,'^:',color='#888',lw=1.3,ms=6,alpha=.8,label='own modular-shift gain',zorder=2)
    ax.plot(x,hed,'v:',color='#2ca25f',lw=1.3,ms=6,alpha=.8,label='hedge factor',zorder=2)
    ax.axhline(1,color='gray',ls=':',lw=1)
    ax.set_xticks(x); ax.set_xticklabels([str(qq) for qq in q])
    ax.set_xlabel(r'prime $q$ dividing the left centre $N$',fontsize=11)
    ax.set_ylabel(r'right-centre shield  $S(d,q)$',fontsize=10.5)
    ax.set_title(f'$6\\Delta N={gap}$  (max error {errlbl})',fontsize=12)
    ax.legend(fontsize=8.5,loc='best'); ax.grid(alpha=.25)
plt.suptitle('Two-centre shield law vs. measurement in $S_{10}$ (23,988,173 twin centres)',fontsize=12.5,y=1.02)
plt.tight_layout()
plt.savefig('fig_paper5_twocentre.pdf',bbox_inches='tight')
plt.savefig('fig_paper5_twocentre.png',dpi=160,bbox_inches='tight')
print("figure saved")
