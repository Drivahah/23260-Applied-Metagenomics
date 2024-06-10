import os 
import pandas as pd
from pathlib import Path
import pycoda

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

#Calculating compositional parameters (resf)
closure_resf = abundance_tbl.coda.closure(1)
#FPKM_resf = 
#log_FPKM_resf =
CLR_resf = abundance_tbl.coda.clr()
#center_resf = 
#variance_resf = abundance_tbl.coda.totvar()




