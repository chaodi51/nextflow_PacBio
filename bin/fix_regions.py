#!/usr/bin/env python3
import sys
import pandas as pd

chrom_dict = {
'pFC17_ss.hSyn1-EF1aint.miR6433.RBG':'pFC17_ss_hSyn1-EF1aint_miR6433_RBG',
'pCC18_ss.CAG.mmc11.RBG':'pCC18_ss_CAG_mmc11_RBG',
'pCC16_ss.CAG.miR6433.RBG':'pCC16_ss_CAG_miR6433_RBG',
'Spk-10001 (BLG-Aldevron)':'Spk-10001-BLG-Aldevron',
'Spk-10001 (BLG-Aldevron-SelfComp)':'Spk-10001-BLG-Aldevron-SelfComp',
'pCC16.ss.CAG.miR6433.RBG':'pCC16_ss_CAG_miR6433_RBG',
'Spark100PK (true)':'Spark100PK',
'Spk-10001-SelfAnnealingModel':'Spk-10001-SelfAnnealingModel',
'Spk-10001-BLG-Aldevron-SelfCompV2':'Spk-10001-BLG-Aldevron-SelfCompV2',
'CNPase.Cre.P2A.EGFP.v2 BGH228 VECTOR':'CNPase-Cre-P2A-EGFP-v2-BGH228',
'MBP.Cre.P2A.EGFP.v2 BGH228 VECTOR':'MBP-Cre-P2A-EGFP-v2-BGH228',
'pCC20_ss.CAG.miR5155':'pCC20-ss-CAG-miR5155',
'pCC41_EF1a.155.miR5155':'pCC41-EF1a-155-miR5155',
'pCC42_EF1a.26.miR5155':'pCC42-EF1a-26-miR5155',
'pCC43_EF1a.33.miR5155':'pCC43-EF1a-33-miR5155',
'VECTOR - pAAV2-CAG_Firefly-Luciferase':'pAAV2-CAG_Firefly-Luciferase'}

region_dict = {
    'Kanamycin resistance':'KanR',
    'Kanamycin-resistance':'KanR',
    'Kan R':'KanR',
    'kan2 marker':'KanR',
    'Kan-R': 'KanR',
    'Spark 100':'Spark100',
    'AAV2 Rep':'AAV2Rep',
    'ApoE-Enhancer-(217-370)': "ApoE-Enhancer"
}

def fix_region(row):
    c = row['region']
    if c in region_dict:
        c = region_dict[c]
    return c

def fix_chrom(row):
    c = row['chrom']
    if c in chrom_dict:
        c = chrom_dict[c]
    return c

in_ls = sys.argv[1:-1]
print(in_ls)
out = sys.argv[-1]
dfs = [pd.read_csv(f, sep='\t') for f in in_ls]
df = pd.concat(dfs)
print(set(df['chrom']))
print( 'debug', len(df[df.chrom=='Spk-10001-SelfAnnealingModel']))
df_filtered1 = df[df['end']-df['st']>20]
print( 'debug1', len(df_filtered1[df_filtered1.chrom=='Spk-10001-SelfAnnealingModel']))
df_filtered2 = df_filtered1[~df_filtered1['region'].str.match('M2s|FNF',na=False)]
print( 'debug2', len(df_filtered2[df_filtered2.chrom=='Spk-10001-SelfAnnealingModel']))
df_filtered2.loc[:, 'chrom'] = df_filtered2.apply(fix_chrom, axis=1)
df_filtered2.loc[:, 'region'] = df_filtered2.apply(fix_region, axis=1)

## simplify regions
contain_yes = "GLA-IDT-co4|HBB2m1|SERPING1co21|ApoE-HCR-1|AAV2|AdV|Ad5|ITR|Promoter|PA|Kan" 
contain_no = "FIX39_Kan|BGH228-PolyA"

df_filtered2 = df_filtered2[df_filtered2.region.str.contains(contain_yes, na=False)]
df_filtered2 = df_filtered2[~df_filtered2.region.str.contains(contain_no, na=False)]
df_filtered2 = df_filtered2.drop_duplicates()
df_filtered2.to_csv(out, sep='\t', index=False)

