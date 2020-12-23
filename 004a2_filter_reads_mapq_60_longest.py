import numpy as np
import pandas as pd
import sys

MAPQ_CUTOFF = 50

col_names_minimap2 = [ str( e ) for e in range( 1, 19 ) ]
mg_num = int( sys.argv[ 1 ] )
print( 'NOW PROCESSING D19-260' + str(mg_num) + '. . .' )

df = pd.read_csv( 'D19-260' + str(mg_num) + '_minimap2_axsr_results.paf', sep='\t', header=None, names=col_names_minimap2 )

d = df.loc[df['12'] >= MAPQ_CUTOFF]
a = d['1'].values

u, c = np.unique(a, return_counts=True)
dup = u[c > 1]

filt_d = d.copy()

repeats = df.loc[df['1'].isin(dup)]
filt_d = pd.concat([filt_d, repeats, repeats]).drop_duplicates(keep=False)

to_add = pd.DataFrame()
for val in dup:
    entries = d.iloc[np.where(a == val)[0]]
    discards = entries.iloc[np.where(entries['12'].values == max(entries[ '12' ].values))[0]]
    discards2 = entries.iloc[np.where(discards['10'].values == max(discards[ '10' ].values))[0]]
    to_add = pd.concat( [ to_add, discards2 ] )

filt_d = pd.concat([filt_d, to_add])

filt_d.to_csv( 'D19-260' + str(mg_num) + '_minimap2_axsr_results_filt.paf', sep='\t', index=False, header=False )
print( 'FILTERING D19-260' + str(mg_num) + ' DONE' )
