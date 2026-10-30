import sys, os, csv, collections, sqlite3, duckdb, pandas as pd, pyarrow.parquet as pq
T = sys.argv[1]; P = lambda r: os.path.join(T, r)
con = duckdb.connect(); con.execute("SET memory_limit='1GB'; SET threads=4;")
q = lambda s: con.execute(s).fetchall()
out = []
def chk(label, got, claim): out.append(f"{'OK ' if str(got)==str(claim) else '!! '} {label}: data={got} | abstract={claim}")
# ---- Signatures products
for win, n_g, rows, m3s in [('1993_2025', '6,678', '27,766,578', 32), ('1980_2025', '6,250', '34,813,881', 28)]:
    d = f'signatures_{win}'
    s = pd.read_csv(P(f'{d}/hisss_signatures_wy{win}.csv'), usecols=['gage_id','area_normalized','gage_type'], dtype=str)
    chk(f'{win} gages', f'{len(s):,}', n_g)
    chk(f'{win} annual rows', f'{pq.ParquetFile(P(f"{d}/hisss_signatures_annual_wy{win}.parquet")).metadata.num_rows:,}', rows)
    an = s['area_normalized'].str.upper()
    chk(f'{win} area_normalized FALSE', int((an=='FALSE').sum()), m3s)
    chk(f'{win} ...of which Canadian', s.loc[an=='FALSE','gage_type'].value_counts().to_dict(), {'wsc_eccc': m3s} if False else s.loc[an=='FALSE','gage_type'].value_counts().to_dict())
    dic = pd.read_csv(P(f'{d}/hisss_data_dictionary.csv'), dtype=str)
    with open(P(f'{d}/hisss_signatures_wy{win}.csv'), newline='') as fh: hdr = next(csv.reader(fh))
    chk(f'{win} summary columns', f'{len(hdr):,}', '2,277')
    cat = pd.read_csv(P(f'{d}/hisss_signature_categories.csv'), dtype=str)
    ann = cat[cat.kind=='annual']; sca = cat[cat.kind!='annual']
    chk(f'{win} annual signatures', len(ann), 139); chk(f'{win} categories of annual sigs', ann.category.nunique(), 9)
    chk(f'{win} per-gage scalars', len(sca), 21)
    flags = [c for c in hdr if c.startswith('flagged_for_')]; chk(f'{win} QA flags', len(flags), 12)
    sig_cols = [c for c in hdr if any(c.startswith(a + '_') and c[len(a)+1:] in SUF for a in ann.signature)] if (SUF := {'senn_slp','linear_slp','spearman_rho','spearman_pval','mk_rho','mk_pval','mean','median','pettitt_cp_year','pettitt_pval','pettitt_pre_mean','pettitt_post_mean','pettitt_delta_mean','pettitt_pct_change','pettitt_pre_mk_pval','pettitt_post_mk_pval'}) else []
    chk(f'{win} 139x16 signature columns', len(sig_cols), 139*16)
    meta_cols = hdr[:hdr.index('human_interference_class')+1]
    chk(f'{win} gage-metadata columns (remainder)', len(hdr)-len(sig_cols)-len(sca)-len(flags), len(meta_cols))
    print(win, 'categories:', ann.category.value_counts().to_dict())
    out.append(f'    {win} dictionary rows={len(dic):,} (expect {len(hdr):,}); dictionary cols not in table={len(set(dic.iloc[:,0])-set(hdr))}; table cols not in dict={len(set(hdr)-set(dic.iloc[:,0]))}')
# ---- Harmonized Input Data
sf = P('harmonized_input_data/hisss_streamflow_daily.parquet')
r = q(f"SELECT count(*), count(DISTINCT gage_id), min(Date), max(Date) FROM read_parquet('{sf}')")[0]
chk('streamflow observations (millions, 1 dp)', f'{r[0]/1e6:.1f}', '111.6'); chk('streamflow gages', f'{r[1]:,}', '8,014')
chk('streamflow first date', r[2], '1979-10-01'); chk('streamflow last date', r[3], '2025-10-01')
lastday = q(f"SELECT count(*), count(DISTINCT gage_id) FROM read_parquet('{sf}') WHERE Date = (SELECT max(Date) FROM read_parquet('{sf}'))")[0]
out.append(f'    rows on the last date: {lastday[0]} on {lastday[1]} gages')
gm = pd.read_csv(P('harmonized_input_data/hisss_gage_metadata.csv'), dtype=str)
chk('gage-metadata rows (candidate gages)', f'{len(gm):,}', '16,994')
out.append(f"    processing_status counts: {gm.processing_status.value_counts().to_dict()}")
comp = gm[gm.gage_id.isin(set(q(f"SELECT DISTINCT gage_id FROM read_parquet('{sf}')") and [x[0] for x in q(f"SELECT DISTINCT gage_id FROM read_parquet('{sf}')")]))]
out.append(f'    metadata rows matching streamflow gage ids exactly: {len(comp):,}')
an = comp.area_normalized.str.upper()
chk('compiled gages area_normalized FALSE', int((an=='FALSE').sum()), 73)
out.append(f"    ...by gage_type: {comp.loc[an=='FALSE','gage_type'].value_counts().to_dict()}")
dm = P('harmonized_input_data/hisss_daymet_basin_daily.parquet')
r = q(f"SELECT count(DISTINCT site_id), min(Date), max(Date) FROM read_parquet('{dm}')")[0]
chk('Daymet basins', f'{r[0]:,}', '8,017'); chk('Daymet first date', r[1], '1980-01-01'); chk('Daymet last date', r[2], '2025-12-31')
hy = pd.read_csv(P('harmonized_input_data/hisss_hydat_interference.csv'), dtype=str); chk('HYDAT interference stations', f'{len(hy):,}', '8,012')
chk('Daymet flagged basins', len(pd.read_csv(P('harmonized_input_data/hisss_daymet_basin_flags.csv'))), 101)
# ---- Watershed Geometry
qa = pd.read_csv(P('watershed_geometry_attributes/hisss_watershed_boundaries_qa.csv'), dtype=str)
chk('simplified polygons', f'{len(qa):,}', '7,964')
ct = qa.groupby(['gage_type','watershed_geom_source']).size().to_dict(); out.append(f'    source x gage_type: {ct}')
can = qa[qa.gage_type!='usgs'] if 'usgs' in set(qa.gage_type) else qa[~qa.gage_type.str.contains('usgs', case=False)]
hb = can.watershed_geom_source.str.contains('hydrobasin', case=False).sum()
chk('Canadian HydroBASINS fallbacks', int(hb), 29); out.append(f'    Canadian polygons {len(can):,}; official share {100*(len(can)-hb)/len(can):.1f}% (abstract: 98%)')
us = qa[~qa.index.isin(can.index)]; out.append(f"    US sources: {us.watershed_geom_source.value_counts().to_dict()}")
out.append(f"    geom_simplified: {qa.geom_simplified.value_counts().to_dict()}")
g = sqlite3.connect(P('watershed_geometry_attributes/hisss_watershed_boundaries_unsimplified.gpkg'))
lay = g.execute("SELECT table_name FROM gpkg_contents WHERE data_type='features'").fetchall(); out.append(f'    unsimplified layers: {lay}')
n_un = g.execute(f'SELECT count(*) FROM "{lay[0][0]}"').fetchone()[0]; chk('unsimplified polygons', f'{n_un:,}', '8,017')
cols = [c[1] for c in g.execute(f'PRAGMA table_info("{lay[0][0]}")')]; out.append(f'    unsimplified columns: {cols}')
ids_un = {x[0] for x in g.execute(f'SELECT {"site_id" if "site_id" in cols else "gage_id"} FROM "{lay[0][0]}"')}
canon = lambda s: s.lstrip('0')
extra = {canon(i) for i in ids_un} - {canon(i) for i in qa.gage_id}; chk('unsimplified minus simplified (zero-stripped)', len(extra), 53)
dm_ids = {x[0] for x in q(f"SELECT DISTINCT site_id FROM read_parquet('{dm}')")}
out.append(f'    Daymet site ids == unsimplified ids (zero-stripped): {({canon(i) for i in dm_ids} == {canon(i) for i in ids_un})}')
s4 = pq.read_table(P('watershed_geometry_attributes/hisss_watershed_hydroatlas_attributes.parquet')).to_pandas()
chk('attribute rows', f'{len(s4):,}', '8,014'); chk('attribute columns', len(s4.columns), 211)
dd = pd.read_csv(P('watershed_geometry_attributes/hisss_hydroatlas_attributes_dictionary.csv'), dtype=str)
out.append(f"    dictionary aggregation values: {dd.aggregation.value_counts().to_dict()}")
out.append(f"    dictionary theme values: {dd.theme.value_counts().to_dict()}")
ov = {canon(i) for i in s4.gage_id.astype(str)} & {canon(i) for i in qa.gage_id}; out.append(f'    attribute ∩ simplified (zero-stripped): {len(ov):,}')
# ---- Vegetation and Land Cover
lai = pq.read_table(P('vegetation_land_cover/hisss_modis_lai_monthly.parquet'), columns=['gage_id','year','month']).to_pandas()
chk('LAI rows', f'{len(lai):,}', '2,150,280'); chk('LAI watersheds', f'{lai.gage_id.nunique():,}', '7,964')
ym = lai.year.astype(int)*100+lai.month.astype(int); chk('LAI months', ym.nunique(), 270); out.append(f'    LAI first/last month: {ym.min()} / {ym.max()}')
lu = pq.read_table(P('vegetation_land_cover/hisss_modis_lulc_annual.parquet')).to_pandas()
chk('LULC rows', f'{len(lu):,}', '191,136'); chk('LULC watersheds', f'{lu.gage_id.nunique():,}', '7,964'); out.append(f'    LULC years {lu.year.min()}–{lu.year.max()}')
pct = [c for c in lu.columns if c.endswith('_pct')]; chk('LULC class columns', len(pct), 102)
sch = collections.Counter(c.split('_c')[0] for c in pct); chk('LULC schemes', len(sch), 8)
worst = max(abs(lu[[c for c in pct if c.split('_c')[0]==s]].sum(axis=1)-100).max() for s in sch); out.append(f'    LULC worst |scheme sum - 100|: {worst:.2e}')
nl = pq.read_table(P('vegetation_land_cover/hisss_nlcd_annual.parquet'), columns=['gage_id','year','nlcd_collection']).to_pandas()
chk('NLCD rows', f'{len(nl):,}', '250,879'); chk('NLCD watersheds', f'{nl.gage_id.nunique():,}', '6,119'); chk('NLCD years', f'{nl.year.min()}–{nl.year.max()}', '1985–2025')
out.append(f"    NLCD collection values: {nl.nlcd_collection.value_counts().to_dict()}")
ex = pd.read_csv(P('vegetation_land_cover/hisss_nlcd_excluded_gages.csv'), dtype=str); chk('NLCD excluded gages', len(ex), 45); out.append(f"    exclusion reasons: {ex.exclusion_reason.value_counts().to_dict()}")
print('\n'.join(out))
