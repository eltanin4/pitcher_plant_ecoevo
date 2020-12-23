import pandas as pd
import numpy as np
import sys

mg_num = int( sys.argv[1] )
bg = pd.read_csv( '~/lte/metagenomes_lte/bedgraphs/bgraph_cov_mg_' + 
                  (3 - len(str(mg_num))) * '0' + str(mg_num) + '.tsv', 
                  sep='\t', header=None, index_col=None )

# Settting parameters and initialising saved quantitites.
avg_contig = {}
unf_contig = {}
rra_genome = {}
chunk_size = 1000   # Relevant for measuring uniformity of coverage.
chunk_cov = []

# Setting up genome names.
g_num_list = list( range( 1, 25 ) ) + list( range( 26, 35 ) )

# Running loop over all genomes.
for g_num in g_num_list:
    g_arr = np.array([int(e[8:10]) for e in bg[0].values]) 

    gf = bg.iloc[ np.where( g_arr == g_num )[0] ]
    rra_genome[ g_num ] = len( gf ) / len( bg )

    len_arr = [int(e.split('_')[5]) for e in gf[0].values]
    gf['len'] = len_arr

    node_arr = [int(e.split('_')[3]) for e in gf[0].values]
    gf['node'] = node_arr

    cov_contig = []
    lens = []
    nodes = np.array(list(set(gf['node'].values)))
    for cn in nodes:
        clen = gf.loc[ gf['node'] == cn ].iloc[0]['len']
        # Filtering out contigs smaller than 1000 bp.
        if clen >= 1000:
            rdf = gf.loc[ gf['node'] == cn ]
            cov_contig.append( sum( ( rdf[2] - rdf[1] ) * rdf[3] ) )
            lens.append( clen )
            num_chunks = int( clen / chunk_size )

            # Measuring average coverage per chunk of size chunk_size. 
            for idxc in range( num_chunks ):
                c_str_pos = chunk_size * idxc 
                c_end_pos = chunk_size * ( idxc + 1 )
                chunk_rdf = rdf.loc[ rdf[2] >= c_str_pos ]
                chunk_rdf = chunk_rdf.loc[ chunk_rdf[2] <= c_end_pos ]
                chunk_cov.append( sum( ( chunk_rdf[2] - chunk_rdf[1] ) * chunk_rdf[3] ) / chunk_size )

    cov_contig = np.array( cov_contig )
    avg_contig[ g_num ] = np.mean( sum( cov_contig ) / sum( lens ) )
    unf_contig[ g_num ] = np.std( chunk_cov )

pd.DataFrame( [ rra_genome, avg_contig, unf_contig ] ).transpose().to_csv( 
                '~/lte/metagenomes_lte/av_covs/mg_' + 
                (3 - len(str(mg_num))) * '0' + str(mg_num) + '_avcov.tsv' , 
                sep='\t', index=False, header=None )

