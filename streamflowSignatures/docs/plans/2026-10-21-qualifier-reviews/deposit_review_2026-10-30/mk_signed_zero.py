import sys, os, math, numpy as np, pandas as pd, pyarrow.parquet as pq
T = sys.argv[1]
def taus(y):
    n = len(y); S = 0
    for i in range(n):
        d = y[i+1:] - y[i]; S += int((d > 0).sum() - (d < 0).sum())
    npairs = n*(n-1)/2
    # correct: one tie group per numeric value (-0.0 == 0.0)
    _, cnt = np.unique(np.where(y == 0, 0.0, y), return_counts=True)
    Ty = sum(t*(t-1)/2 for t in cnt)
    # buggy (Julia unique/isequal + count(==)): zero group counted once per distinct signed zero
    nz = int((y == 0).sum()); n_signed = int(np.any((y == 0) & ~np.signbit(y))) + int(np.any((y == 0) & np.signbit(y)))
    Ty_bug = Ty + (n_signed - 1) * nz*(nz-1)/2 if nz > 1 else Ty
    f = lambda T_: S / math.sqrt((npairs - T_) * npairs) if npairs - T_ > 0 else float('nan')
    return f(Ty), f(Ty_bug)
for win in ['1993_2025', '1980_2025']:
    pf = pq.ParquetFile(os.path.join(T, f'signatures_{win}/hisss_signatures_annual_wy{win}.parquet'))
    mixed = []
    for rg in range(pf.num_row_groups):
        t = pf.read_row_group(rg).to_pandas()
        z = t[t.value == 0]
        if z.empty: continue
        g = z.assign(neg=np.signbit(z.value.to_numpy())).groupby(['gage_id','signature']).neg.agg(['min','max'])
        keys = g[(g['min'] == False) & (g['max'] == True)].index
        if len(keys):
            sub = t.set_index(['gage_id','signature']).loc[keys].reset_index()
            mixed.append(sub)
    mixed = pd.concat(mixed) if mixed else pd.DataFrame(columns=['gage_id','signature','water_year','value'])
    pairs = mixed[['gage_id','signature']].drop_duplicates()
    print(f'{win}: series holding both +0.0 and -0.0: {len(pairs)} on {pairs.gage_id.nunique()} gages; signatures: {pairs.signature.value_counts().to_dict()}')
    if pairs.empty: continue
    cols = ['gage_id'] + sorted({f'{s}_mk_rho' for s in pairs.signature})
    summ = pd.read_csv(os.path.join(T, f'signatures_{win}/hisss_signatures_wy{win}.csv'), usecols=cols, dtype={'gage_id': str}).set_index('gage_id')
    n_cmp = n_bug = n_ok = 0
    for (gid, sig), grp in mixed.groupby(['gage_id','signature']):
        y = grp.sort_values('water_year').value.to_numpy(); y = y[np.isfinite(y)]
        dep = summ.at[gid, f'{sig}_mk_rho']
        if not np.isfinite(dep): continue
        good, bug = taus(y); n_cmp += 1
        if abs(dep - bug) < 1e-9 and abs(dep - good) > 1e-9: n_bug += 1
        elif abs(dep - good) < 1e-9: n_ok += 1
        if n_cmp <= 2: print(f'   e.g. {gid} {sig}: deposited {dep:.10f}  buggy-ties {bug:.10f}  correct {good:.10f}')
    print(f'   deposited mk_rho compared: {n_cmp}; reproduce the double-counted ties: {n_bug}; match correct ties: {n_ok}')
