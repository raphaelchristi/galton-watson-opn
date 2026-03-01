#!/usr/bin/env python3
"""
Galton-Watson Analysis of Cyclotomic Cascades (FAST)
====================================================
Optimized: ~10-20s on Mac Mini M4.

Key optimizations vs v1:
  - 300 seeds (not 850)
  - Trial division only up to 10K (skip Pollard rho entirely)
  - Max 8 rounds, max 200 primes per cascade
  - Cached Φ₃ factorizations
  - No scatter plot data collection

Author: Raphael Christi
Date: March 2026
"""

import math, time
from functools import lru_cache
from collections import defaultdict

# ============================================================
# PRIMES
# ============================================================

def _sieve(n):
    s = bytearray(b'\x01') * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i*i::i] = b'\x00' * len(s[i*i::i])
    return [i for i, v in enumerate(s) if v]

PRIMES = _sieve(100000)
P_SET = set(PRIMES)

def is_prime(n):
    if n <= 100000: return n in P_SET
    if n < 2: return False
    if any(n % p == 0 for p in PRIMES if p * p <= n and p <= 1000):
        return False
    d, r = n - 1, 0
    while d % 2 == 0: d //= 2; r += 1
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if a >= n: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

# ============================================================
# Φ₃ FACTORIZATION (cached)
# ============================================================

@lru_cache(maxsize=50000)
def phi3_factors(p):
    """Odd prime factors of Φ₃(p) = p²+p+1. Trial div to 10K + primality check on cofactor."""
    n = p * p + p + 1
    while n % 2 == 0: n //= 2
    factors = []
    for q in PRIMES:
        if q > 10000 or q * q > n: break
        if n % q == 0:
            factors.append(q)
            while n % q == 0: n //= q
    if n > 1 and is_prime(n):
        factors.append(n)
    # if n > 1 and not prime: composite cofactor, discard (rare, large)
    return tuple(factors)

# ============================================================
# OFFSPRING COMPUTATION
# ============================================================

def compute_all_offspring(seeds, max_rounds=8, max_size=200):
    """
    Run cascades and collect (p, n_children) for every prime encountered.
    Returns: list of (p, n_children), n_closed, n_diverged
    """
    records = []
    n_closed = 0
    n_diverged = 0

    for seed in seeds:
        S = set(seed)
        frontier = sorted(seed)
        closed = False

        for _ in range(max_rounds):
            next_frontier = []
            for p in frontier:
                children = [q for q in phi3_factors(p) if q not in S and q > 2]
                records.append((p, len(children)))
                S.update(children)
                next_frontier.extend(children)

            if not next_frontier:
                closed = True
                break
            frontier = sorted(set(next_frontier))
            if len(S) > max_size:
                break

        if closed:
            n_closed += 1
        else:
            n_diverged += 1

    return records, n_closed, n_diverged

# ============================================================
# MAIN
# ============================================================

def main():
    t0 = time.time()

    print("=" * 70)
    print("GALTON-WATSON ANALYSIS OF CYCLOTOMIC CASCADES")
    print("=" * 70)

    # ── Seeds ───────────────────────────────────────────────
    pm3 = [p for p in PRIMES if p > 3 and p % 3 == 1 and p < 10000]

    seeds = (
        [(p,) for p in pm3[:200]] +                       # 200 singles
        [(3, p) for p in pm3[:60]] +                       # 60 pairs
        [(3, p, q) for i, p in enumerate(pm3[:20])
                    for q in pm3[i+1:i+3]]                 # 40 triples
    )
    print(f"\n  Seeds: {len(seeds)}  (singles + pairs + triples)")

    # ── Run ─────────────────────────────────────────────────
    records, n_closed, n_diverged = compute_all_offspring(seeds)
    t_run = time.time() - t0
    print(f"  Cascades done in {t_run:.1f}s — {len(records):,} offspring records")
    print(f"  Closed: {n_closed}   Diverged: {n_diverged}")

    # ── Global μ ────────────────────────────────────────────
    children = [r[1] for r in records]
    mu = sum(children) / len(children)

    print(f"\n{'─' * 70}")
    print(f"  GLOBAL  μ̄ = {mu:.4f}  ", end="")
    if mu > 1:
        print("★ SUPERCRITICAL")
    else:
        print("subcritical")
    print(f"{'─' * 70}")

    # ── Distribution ────────────────────────────────────────
    dist = defaultdict(int)
    for c in children: dist[c] += 1

    print(f"\n  Offspring distribution:")
    print(f"  {'k':>4} {'count':>8} {'frac':>8}  ")
    for k in sorted(dist)[:13]:
        f = dist[k] / len(children)
        print(f"  {k:>4} {dist[k]:>8} {f:>8.3f}  {'█' * int(f * 60)}")
    rest = sum(v for k, v in dist.items() if k > 12)
    if rest:
        print(f"  {'13+':>4} {rest:>8} {rest/len(children):>8.3f}")

    # ── μ(p) by log-bin ────────────────────────────────────
    bins = defaultdict(list)
    for p, nc in records:
        b = round(math.log10(max(p, 2)) * 2) / 2
        bins[b].append(nc)

    print(f"\n  μ(p) by prime size:")
    print(f"  {'log₁₀(p)':>9} {'~range':>14} {'n':>7} {'μ':>7} {'σ':>7} {'sup?':>5}")
    print(f"  {'─'*9} {'─'*14} {'─'*7} {'─'*7} {'─'*7} {'─'*5}")

    bin_data = []
    for b in sorted(bins):
        v = bins[b]
        n = len(v)
        m = sum(v) / n
        s = (sum((x - m)**2 for x in v) / max(n - 1, 1)) ** 0.5
        lo, hi = int(10**(b - .25)), int(10**(b + .25))
        sup = "★" if m > 1 else ""
        print(f"  {b:>9.1f} {lo:>6}-{hi:<6} {n:>7} {m:>7.3f} {s:>7.3f} {sup:>5}")
        bin_data.append((b, m, s, n))

    # ── μ by residue class ──────────────────────────────────
    by_mod = defaultdict(list)
    for p, nc in records: by_mod[p % 3].append(nc)

    print(f"\n  μ by residue (mod 3):")
    for c in sorted(by_mod):
        v = by_mod[c]
        m = sum(v) / len(v)
        print(f"    p ≡ {c} (mod 3): μ = {m:.3f}  (n={len(v)})")

    # ── Largest child / parent ──────────────────────────────
    ratios = []
    for p in pm3[:300]:
        facs = phi3_factors(p)
        if facs:
            biggest = max(facs)
            if p > 5:
                ratios.append(math.log10(biggest) / math.log10(p))

    if ratios:
        med = sorted(ratios)[len(ratios)//2]
        avg = sum(ratios) / len(ratios)
        print(f"\n  Largest child digit ratio (log(max_child)/log(p)):")
        print(f"    Mean: {avg:.3f}   Median: {med:.3f}   (Stewart → 2.0)")

    # ── Multi-Φ comparison ──────────────────────────────────
    print(f"\n  μ for different cyclotomic polynomials:")
    for name, deg, fn in [
        ("Φ₃", 2, lambda p: p*p + p + 1),
        ("Φ₅", 4, lambda p: p**4 + p**3 + p**2 + p + 1),
        ("Φ₇", 6, lambda p: p**6 + p**5 + p**4 + p**3 + p**2 + p + 1),
    ]:
        tc = tp = 0
        for p in pm3[:150]:
            val = fn(p)
            while val % 2 == 0: val //= 2
            nf = 0
            tmp = val
            for q in PRIMES:
                if q > 10000 or q * q > tmp: break
                if tmp % q == 0:
                    nf += 1
                    while tmp % q == 0: tmp //= q
            if tmp > 1 and is_prime(tmp): nf += 1
            tc += nf; tp += 1
        print(f"    {name} (deg {deg}): μ = {tc/tp:.3f}")

    # ── Save HTML plot ──────────────────────────────────────
    labels = [f"{b:.1f}" for b, _, _, _ in bin_data]
    mus = [round(m, 4) for _, m, _, _ in bin_data]
    dl = [str(k) for k in sorted(dist) if k <= 12]
    dv = [round(dist[k]/len(children), 4) for k in sorted(dist) if k <= 12]

    html = f"""<!DOCTYPE html>
<html><head><title>Galton-Watson — Φ₃ Cascades</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0"></script>
<style>
body{{font-family:system-ui;max-width:1100px;margin:40px auto;padding:0 20px;background:#0a0a0a;color:#e0e0e0}}
h1{{color:#4fc3f7}}h2{{color:#81c784}}.box{{background:#1a1a1a;border-radius:12px;padding:20px;margin:20px 0}}
canvas{{max-height:380px}}.result{{font-size:2em;color:#ff9800;text-align:center;padding:20px;background:#1a1a1a;border-radius:12px;border:2px solid #ff9800;margin:20px 0}}
.g{{display:grid;grid-template-columns:1fr 1fr;gap:20px}}@media(max-width:768px){{.g{{grid-template-columns:1fr}}}}
</style></head><body>
<h1>Galton-Watson Analysis — Cyclotomic Cascades</h1>
<div class="result">μ̄ = {mu:.4f} — {'★ SUPERCRITICAL' if mu > 1 else 'SUBCRITICAL'}</div>
<div class="g">
<div class="box"><h2>μ(p) vs log₁₀(p)</h2><canvas id="c1"></canvas></div>
<div class="box"><h2>Offspring Distribution</h2><canvas id="c2"></canvas></div>
</div>
<script>
new Chart(document.getElementById('c1'),{{type:'bar',data:{{labels:{labels},datasets:[
{{label:'μ(p)',data:{mus},backgroundColor:'rgba(79,195,247,0.7)',borderColor:'#4fc3f7',borderWidth:1}},
{{label:'μ=1 critical',data:Array({len(labels)}).fill(1),type:'line',borderColor:'#ff5252',borderDash:[5,5],pointRadius:0,borderWidth:2,fill:false}}
]}},options:{{scales:{{x:{{title:{{display:true,text:'log₁₀(p)',color:'#aaa'}},ticks:{{color:'#aaa'}},grid:{{color:'#333'}}}},y:{{title:{{display:true,text:'μ(p)',color:'#aaa'}},ticks:{{color:'#aaa'}},grid:{{color:'#333'}},min:0}}}},plugins:{{legend:{{labels:{{color:'#ccc'}}}}}}}}}});
new Chart(document.getElementById('c2'),{{type:'bar',data:{{labels:{dl},datasets:[{{label:'Fraction',data:{dv},backgroundColor:'rgba(129,199,132,0.7)',borderColor:'#81c784',borderWidth:1}}]}},options:{{scales:{{x:{{title:{{display:true,text:'Children',color:'#aaa'}},ticks:{{color:'#aaa'}},grid:{{color:'#333'}}}},y:{{title:{{display:true,text:'Fraction',color:'#aaa'}},ticks:{{color:'#aaa'}},grid:{{color:'#333'}}}}}},plugins:{{legend:{{labels:{{color:'#ccc'}}}}}}}}}});
</script></body></html>"""

    with open("galton_watson_plot.html", "w") as f:
        f.write(html)

    # ── Final ───────────────────────────────────────────────
    elapsed = time.time() - t0
    print(f"\n{'=' * 70}")
    print(f"  μ̄ = {mu:.4f}  {'★ SUPERCRITICAL' if mu > 1 else '⚠ subcritical'}")
    print(f"  Runtime: {elapsed:.1f}s")
    print(f"  Output: galton_watson_plot.html")
    print(f"{'=' * 70}")

if __name__ == "__main__":
    main()
