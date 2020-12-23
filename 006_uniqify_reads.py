import numpy as np
import pandas as pd
import sys

# Reading the file from the arugment.
col_names_minimap2 = [ str( e ) for e in range( 1, 19 ) ]
mg_num = int( sys.argv[ 1 ] )
print( 'NOW PROCESSING D19-2600' + str(mg_num) + '. . .' )

df = pd.read_csv( 'D19-2600' + str(mg_num) + '_minimap2_axsr_results_filt.paf', sep='\t', header=None, names=col_names_minimap2 )

# Getting the read names and removing the forward-reverse information.
a = df[ '1' ].values
a = np.array([ e[:-2] for e in a ])

# Figuring out duplicates.
u, c = np.unique(a, return_counts=True)
dup = u[c > 2]

# Removing reverse /2 reads for the duplicates.
filt_d = df.copy()
rems = []
for t in dup:
    rems.append( t + '/2' )

filt_d = filt_d[ ~filt_d[ '1' ].isin(rems) ]

# Finding out reads which still map to more than contig with the same quality and length.
u, c = np.unique(filt_d['1'].values, return_counts=True)
dup2 = u[c > 1]

# Choosing the one that maps to the smallest contig, so the highest alignment percentage.
repeats = filt_d.loc[filt_d['1'].isin(dup2)]
filt_d_final = pd.concat([filt_d, repeats, repeats]).drop_duplicates(keep=False)

to_add = pd.DataFrame()
for val in dup2:
    entries = filt_d.iloc[np.where(filt_d['1'].values == val)[0]]
    discards = entries.iloc[np.where(entries['7'].values == min(entries[ '7' ].values))[0]]
    if len( discards ) > 1:
        discards = discards.sample(1)
    to_add = pd.concat( [ to_add, discards ] )

filt_d_final = pd.concat([filt_d_final, to_add])

# Checking if no more duplicate reads.
u, c = np.unique(filt_d_final['1'].values, return_counts=True)
dup3 = u[c > 1]

# Checking if all reads map to a unique contig now.
if len( dup3 ):
    print( 'COMPLETED, BUT STILL HAVE' + str( len( dup3 ) ) + 'DUPLICATES.')

# Saving new filtered paf file.
filt_d.to_csv( 'D19-2600' + str(mg_num) + '_minimap2_axsr_results_filt_dup.paf', sep='\t', index=False, header=False )
print( 'FILTERING D19-2600' + str(mg_num) + ' DONE' )
