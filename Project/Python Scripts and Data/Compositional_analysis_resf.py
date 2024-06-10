import os 
import pandas as pd
from pathlib import Path
import pycoda
import math

# Path to dir with mapstat files
PATH = os.getcwd()


mapstat_files = [] 
for mfile in Path(PATH).glob('*silva.mapstat'):
    df = pd.read_csv(mfile, skiprows=6, sep='\t')
    cols = list(df.columns)
    # strip trailing hash from columns
    cols[0] = cols[0].replace('# ', '')
    df.columns = cols
    # add sample name to mapstat file
    df['sample_name'] = mfile.name.replace('.mapstat', '')
    mapstat_files.append(df)
# concatinate tables
catmapstat = pd.concat(mapstat_files)
# create a pivot table with samples as rows and features as columns
abundance_tbl = catmapstat.pivot(index='sample_name', columns='refSequence', values='fragmentCountAln')
# replace null values with zeros
abundance_tbl = abundance_tbl.fillna(0)

#Calculating compositional parameters (SILVA)
closure_silva = abundance_tbl.coda.closure(1)
#FPKM_silva = 
#log_FPKM_silva =
CLR_silva = abundance_tbl.coda.clr()
#center_silva = 
#variance_silva = abundance_tbl.coda.totvar()


# Path to dir with mapstat files
PATH = os.getcwd()


mapstat_files = [] 
for mfile in Path(PATH).glob('*_resf.mapstat'):
    df = pd.read_csv(mfile, skiprows=6, sep='\t')
    cols = list(df.columns)
    # strip trailing hash from columns
    cols[0] = cols[0].replace('# ', '')
    df.columns = cols
    # add sample name to mapstat file
    df['sample_name'] = mfile.name.replace('.mapstat', '')
    mapstat_files.append(df)
# concatinate tables
catmapstat = pd.concat(mapstat_files)
# create a pivot table with samples as rows and features as columns
abundance_tbl = catmapstat.pivot(index='sample_name', columns='refSequence', values='fragmentCountAln')
# replace null values with zeros
abundance_tbl = abundance_tbl.fillna(0)
#abundance_tbl = abundance_tbl.coda.zero_replacement(n_samples=5000)

res_gene_length = pd.read_csv('kma_Resfinder.name.length.txt', sep ='\t')

'''
resf_count_len = [] #Dataframes with gene names, read counts and gene lengths

for _ in range(len(mapstat_files)):
    df = mapstat_files[_].merge(res_gene_length, on='refSequence')
    df = df[['refSequence', 'fragmentCountAln', 'Gene_length']]
    df = df.fillna(0)
    df.coda.zero_replacement(n_samples=5000)
    resf_count_len.append(df)
del df
'''

df = abundance_tbl.T.copy()
for _ in df:
    df[_] += 1/(df[_].sum())
resf_count_len = df.merge(res_gene_length, on='refSequence')

del df

M = {'ERR2597508_resf' : 185303, 
     'ERR2597509_resf' : 115735, 
     'ERR2597510_resf' : 163175}

cols = resf_count_len.columns
cols = cols[1:4]

def log_FPKM(mapped_reads, tot_frags, gene_length):
    FPKM = (mapped_reads * 10**9)/(tot_frags * gene_length)
    return math.log(FPKM)

FPKM_table = {'refSequence' : resf_count_len['refSequence'],
              'ERR2597508_resf' : [], 
              'ERR2597509_resf' : [], 
              'ERR2597510_resf' : []}

for _ in range(len(resf_count_len['ERR2597508_resf'])):
    FPKM_table['ERR2597508_resf'].append(log_FPKM(resf_count_len['ERR2597508_resf'][_], M['ERR2597508_resf'], resf_count_len['Gene_length'][_]))  

for _ in range(len(resf_count_len['ERR2597509_resf'])):
    FPKM_table['ERR2597509_resf'].append(log_FPKM(resf_count_len['ERR2597509_resf'][_], M['ERR2597509_resf'], resf_count_len['Gene_length'][_]))  
    
for _ in range(len(resf_count_len['ERR2597510_resf'])):
    FPKM_table['ERR2597510_resf'].append(log_FPKM(resf_count_len['ERR2597510_resf'][_], M['ERR2597510_resf'], resf_count_len['Gene_length'][_]))  
    
FPKM_table = pd.DataFrame.from_dict(FPKM_table)

### CALCULATE SUM OF FPKM FOR EACH SAMPLE

#Calculating compositional parameters (resf)
closure_resf = abundance_tbl.coda.closure(1)
#FPKM_resf = 
#log_FPKM_resf =
CLR_resf = abundance_tbl.coda.clr()
#center_resf = 
#variance_resf = abundance_tbl.coda.totvar()



