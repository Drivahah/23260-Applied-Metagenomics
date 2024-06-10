# %% Libraries
import os 
import pandas as pd
from pathlib import Path
import math
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# %% LOAD MAPSTAT AND OBTAIN ABUNDANCE TABLE

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


# %% LOAD GENE LENGTH FILE AND CALCULATE FPKM
    
res_gene_length = pd.read_csv('kma_Resfinder.name.length.txt', sep ='\t')

df = abundance_tbl.T.copy()
for _ in df:
    df[_] += 1/(df[_].sum())
resf_count_len = df.merge(res_gene_length, on='refSequence')

del df

def log_FPKM(mapped_reads, tot_frags, gene_length):
    FPKM = (mapped_reads * 10**9)/(tot_frags * gene_length)
    return math.log(FPKM, 10)

M = {'ERR2597508_resf' : 185303, #Count of bacteria from SILVA
     'ERR2597509_resf' : 115735, 
     'ERR2597510_resf' : 163175}

FPKM_table = {'refSequence' : resf_count_len['refSequence'],
              'ERR2597508_resf' : [], 
              'ERR2597509_resf' : [], 
              'ERR2597510_resf' : []}

# Load FPKM values in FPKM_table
for _ in range(len(resf_count_len['ERR2597508_resf'])):
    FPKM_table['ERR2597508_resf'].append(log_FPKM(resf_count_len['ERR2597508_resf'][_], M['ERR2597508_resf'], resf_count_len['Gene_length'][_]))  

for _ in range(len(resf_count_len['ERR2597509_resf'])):
    FPKM_table['ERR2597509_resf'].append(log_FPKM(resf_count_len['ERR2597509_resf'][_], M['ERR2597509_resf'], resf_count_len['Gene_length'][_]))  
    
for _ in range(len(resf_count_len['ERR2597510_resf'])):
    FPKM_table['ERR2597510_resf'].append(log_FPKM(resf_count_len['ERR2597510_resf'][_], M['ERR2597510_resf'], resf_count_len['Gene_length'][_]))  
    
FPKM_table = pd.DataFrame.from_dict(FPKM_table)

# Sum of FPKM for each sample
FPKM_total = pd.DataFrame(FPKM_table.loc[:, FPKM_table.columns!='refSequence'].sum())


# %% PLOT TOTAL FPKM
x_pos = np.arange(len(FPKM_total))

# Create bars and choose color
plt.bar(x_pos, FPKM_total[0], color = (0.5,0.1,0.5,0.75))

# Add title and axis names
plt.title('Total log2(FPKM)')
plt.xlabel('Samples')
plt.ylabel('Total log2(FPKM)')

FPKM_total['Rounded'] = FPKM_total[0].round(1)

# Create names on the x axis
plt.xticks(x_pos, list(FPKM_total.index))

# Create labels
label = FPKM_total[0].apply(str)
 
# Text on the top of each bar
for _ in x_pos:
    plt.text(x = _ -0.15, y = FPKM_total[0][_] - 35, s = FPKM_total['Rounded'][_], size = 13, color = 'White')
    
# Save as PDF
plt.savefig('FPKM_total.svg', bbox_inches='tight') 

# Show graph
plt.show()

# %% DISTANCE BETWEEN SAMPLES

distance_08_09 = np.linalg.norm(FPKM_table['ERR2597508_resf'] - FPKM_table['ERR2597509_resf'])
distance_08_10 = np.linalg.norm(FPKM_table['ERR2597508_resf'] - FPKM_table['ERR2597510_resf'])
distance_10_09 = np.linalg.norm(FPKM_table['ERR2597510_resf'] - FPKM_table['ERR2597509_resf'])

# %% PCA

# Properly formatting FPKM_table for PCA
x = FPKM_table.T
gene_names = x.iloc[0]
x = x[1:]
x.columns = gene_names

# Standardizing the features
x = StandardScaler().fit_transform(x)

# Perform PCA
pca = PCA(n_components=2) 
principalComponents = pca.fit_transform(x)
expl_var = pca.explained_variance_ratio_
finalDf = pd.DataFrame(data = principalComponents,
                           columns = ['principal component 1', 'principal component 2'])

'''
I tried also with 3 components, but the 
third one explained almost no variance
'''

# Visualize 2D projection
labels = [f"PC {i+1} ({var:.1f}%)" for i, var in enumerate(expl_var * 100)]
fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel(labels[0], fontsize = 15)
ax.set_ylabel(labels[1], fontsize = 15)
ax.set_title('FPKM PCA', fontsize = 20)
targets = ['08', '09', '10']
finalDf['target'] = targets
colors = ['r', 'g', 'b']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['target'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 100)
ax.legend(targets)
ax.grid()

# Save as PDF
plt.savefig('FPKM_PCA.svg', bbox_inches='tight') 

# Show graph
plt.show()

# %% GROUPING RESISTANCES

res_groups = pd.read_csv('ResFinder_classes.tsv', sep='\t')
FPKM_table = FPKM_table.merge(res_groups, left_on='refSequence', right_on='refSequence')
res_groups_FPKM = FPKM_table.groupby('resGroup').sum()

# %% PLOT GROUPED RESISTANCES
 
# width of the bars
barWidth = 0.3
 
# Choose the height of the bars
bars1 = res_groups_FPKM['ERR2597508_resf']
bars2 = res_groups_FPKM['ERR2597509_resf']
bars3 = res_groups_FPKM['ERR2597510_resf']
 
# The x position of bars
r1 = np.arange(len(bars1))
r2 = [x + barWidth for x in r1]
r3 = [x + 2 * barWidth for x in r1]

# Color mapping
colors = iter([plt.cm.tab10(i) for i in range(10)])

# Create vertical bars
plt.bar(r1, bars1, width = barWidth, color = [next(colors)], alpha = 0.7, capsize=7, label='08')
plt.bar(r2, bars2, width = barWidth, color = [next(colors)], alpha = 0.8, capsize=7, label='09')
plt.bar(r3, bars3, width = barWidth, color = [next(colors)], alpha = 0.9, capsize=7, label='10')

plt.title('log2(FPKM) of resistance genes categories')
plt.xticks([r + barWidth for r in range(len(bars1))], res_groups_FPKM.index, rotation=90)
plt.yticks(ticks=np.arange(0, 200, 25))
plt.ylabel('log2(FPKM)')
plt.legend()

# Save as PDF
plt.savefig('FPKM_resistance_categories.svg', dpi=100, bbox_inches='tight') 

# Show graphic
plt.show()
