# 6N Twin-Prime Two-Centre Shield (Part V)

A measured, quantitative law for the **two-centre correlation** in the twin-prime
gap structure on the 6N ± 1 skeleton: the **modular-shift gain** times a
**cross-prime hedge**.

**Background.** Parts II–IV showed the twin-gap distribution depends on ω₍>3₎(N).
Part IV closed the long-gap (210) residual with a single-centre spatial-
compression penalty but left the short-gap (42) residual, and identified the
genuine **two-centre correlation** between the factorisations of N and N+d as the
open problem. This part measures that correlation directly.

**The object.** For a centre-step d and prime q, the two-centre shield is

```
    S(d,q) = P(N+d twin | q|N, N twin) / P(N+d twin | N twin)
```

measured over a high-ω band (ω₍>3₎ ∈ {4,5,6}) of S₁₀.

**The law (closes to ≤1.5% on 23,988,173 twin centres).**

```
    S(d,q) = g(d,q) · h(d,q)
    g (modular-shift gain): q|N forces N+d ≡ d (mod q); if d is q-safe the right
        centre is certainly q-safe, vs its MEASURED baseline q-safe rate β_q(d).
        g = 1{d ∉ dead(q)} / β_q(d),  dead(q) = {±6⁻¹ mod q}.
    h (cross-prime hedge): conditioning on q|N moves the right centre's residues
        at the OTHER primes q'; h = ∏_{q'≠q} P(N+d q'-safe | q|N) / β_{q'}(d).
```

Both factors require the **measured** modular distribution of twin centres (which
carries the Part-I enrichment); an idealised uniform baseline overstates g and
misses h entirely.

**Two findings that overturn natural guesses.**
- The **42 shield is a q=5 effect, not q=7.** Since N+7 ≡ N+2 (mod 5), having
  5|N places the right centre at the safe residue 2. The prime 7 that divides
  the physical gap 42 is nearly inert (g(42,7)=1).
- The **hedge sign differs by gap**: for 42 it is < 1 (the shift makes the right
  centre collide more often at q=5); for 210 it is > 1. This is why the
  single-centre model of Part IV under-/over-shot.

> **Scope.** Experimental / computational number theory. This is the
> *right-centre conditional survival* law. Bridging it to the ω-stratified gap
> preference r(d|ω) of Part IV — and with it the final fate of the 42 residual —
> is **not done here** and is left as the open problem. No claim is made about
> the infinitude of twin primes or any prime k-tuple conjecture.

Part I: doi:10.5281/zenodo.20470367 · II: doi:10.5281/zenodo.20477664 ·
III: doi:10.5281/zenodo.20498668 · IV: doi:10.5281/zenodo.20500465

---

## Layout

```
.
├── README.md
├── LICENSE                 (MIT)
├── CITATION.cff
├── data/
│   └── tc_plot_S10.csv      gap, q, own_gain, hedge, predicted, measured  (S10)
├── code/
│   ├── interferometer.py        the discovery tool: measures right-centre
│   │                            survival grouped by left-centre divisibility by {5,7}
│   ├── full_twocentre.py        the result: measures S(d,q) and the g×h
│   │                            decomposition for all primes; emits tc_plot_S{K}.csv
│   └── make_twocentre_fig.py    builds the 2-panel prediction-vs-measurement figure
├── figures/                fig_paper5_twocentre.{pdf,png}
└── paper/                  Chen_6N_Paper5.{tex,pdf} + figure
```

## Reproducing

Requirements: Python 3.8+, `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib

# 1. The interferometer (fast discovery view: 42 shield is 5-driven, not 7).
python code/interferometer.py            # S10 (~15-20 min)
MAXK=9 python code/interferometer.py     # S9 (faster)

# 2. The full two-centre law (g × h) for all primes; writes tc_plot_S{K}.csv.
python code/full_twocentre.py            # S10
MAXK=9 python code/full_twocentre.py     # S9

# 3. The figure (reads ../data/tc_plot_S10.csv).
cd code && python make_twocentre_fig.py
```

### Conventions (same as Parts II–IV)

- Twin centre: N with 6N−1, 6N+1 both prime. Centre-step d; physical gap 6d.
- dead(q) = {±6⁻¹ mod q}: a centre is q-safe iff its residue ∉ dead(q).
- Right centre M = N+d; q|N forces M ≡ d (mod q) — the source of the gain.
- High-ω band ω∈{4,5,6} for comparability; primes that lock the gap
  (q | 6d±1: 11,19 for 210; 41,43 for 42) give shield 0 and are excluded from
  the fit table.
- Engine: complete segmented-sieve factorisation + deterministic interval-sieve
  primality; S₁₀ twin count 23,988,173 matches Part I. Membership test of
  "N+d is a twin centre" by binary search on the sorted twin array (memory-light).

## License

MIT — see `LICENSE`.
