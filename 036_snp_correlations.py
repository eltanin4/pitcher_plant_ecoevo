import numpy as np
import pandas as pd
import sys
from scipy.stats import pearsonr, spearmanr
from itertools import combinations
import scipy
import scipy.cluster.hierarchy as sch

m01 = [ 35, 43, 51, 59, 67, 75, 83, 91 ]
m02 = [ 36, 44, 52, 60, 68, 76, 84, 92 ]
m03 = [ 37, 45, 53, 61, 69, 77, 85, 93 ]
m04 = [ 38, 46, 54, 62, 70, 78, 86, 94 ]
m05 = [ 39, 47, 55, 63, 71, 79, 87, 95 ]
m06 = [ 40, 48, 56, 64, 72, 80, 88, 96 ]
m07 = [ 41, 49, 57, 65, 73, 81, 89, 97 ]
m08 = [ 42, 50, 58, 66, 74, 82, 90, 98 ]
m09 = [ 99, 101, 103, 105, 107, 109, 111, 113 ]
m10 = [ 100, 102, 104, 106, 108, 110, 112, 114 ]

allms = [ m01, m02, m03, m04, m05, m06, m07, m08, m09, m10 ]
g_num_list = list( range( 1, 25 ) ) + list( range( 26, 35 ) )

MIN_COV = 20

# this_mic = int( sys.argv[ 1 ] )
max_freqs, fix_times, appear_times = [], [], []
valid_snps, corr_matrix, sign_matrix = {}, {}, {}
for this_mic in range( 1, 11 ):
    tm00 = allms[ this_mic - 1 ]
    print( 'Processing microcosm ' + str( tm00 ) )
    
    # Initializing correlation matrix.
    valid_snps[ this_mic - 1 ] = {}
    corr_matrix[ this_mic - 1 ] = {}
    sign_matrix[ this_mic - 1 ] = {}

    is_invalid = False
    g_num = int( sys.argv[ 1 ] )
    print( 'Processing genome ' + str( g_num ) + '.' )
    try:
        df = pd.read_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/snps/all_alleles/fb_trajs/mic_'
                          + str( this_mic ) + '_gnm_' + str( g_num ) +
                          '_snp_trajectories.tsv', sep='\t', index_col=0 )

        if len( df ) == 0:
            is_invalid = True
            continue
    except:
        is_invalid = True
        continue

    # Gathering valid SNPs.
    valid_snps[ this_mic - 1 ][ g_num ] = []
    for idx, row in df.iterrows():
        # print( idx, row )
        these_depths = [ row[ 't1_locus_depth' ], row[ 't2_locus_depth' ], row[ 't3_locus_depth' ], row[ 't4_locus_depth' ],
                         row[ 't5_locus_depth' ], row[ 't6_locus_depth' ], row[ 't7_locus_depth' ], row[ 't8_locus_depth' ] ]
        these_ref_freqs = [ row[ 't1_ref_freq' ], row[ 't2_ref_freq' ], row[ 't3_ref_freq' ], row[ 't4_ref_freq' ],
                            row[ 't5_ref_freq' ], row[ 't6_ref_freq' ], row[ 't7_ref_freq' ], row[ 't8_ref_freq' ] ]
        these_alt_freqs = [ row[ 't1_alt_freq' ], row[ 't2_alt_freq' ], row[ 't3_alt_freq' ], row[ 't4_alt_freq' ],
                            row[ 't5_alt_freq' ], row[ 't6_alt_freq' ], row[ 't7_alt_freq' ], row[ 't8_alt_freq' ] ]
        if min( these_depths ) > MIN_COV and min( these_alt_freqs ) < 0.95 and these_alt_freqs[ 0 ] < 0.95:
            if len( [ e for e in these_alt_freqs if e not in [ 0.0, 1.0 ] ] ) > 0:
                valid_snps[ this_mic - 1 ][ g_num ].append( these_alt_freqs )

    # Calculating correlation matrix for this genome if more than 50 SNPs for this genome.
    NUM_VALID_SNPS = 20
    if len( valid_snps[ this_mic - 1 ][ g_num ] ) > NUM_VALID_SNPS:
        this_matrix = np.zeros( ( len( valid_snps[ this_mic - 1 ][ g_num ] ), 
                                len( valid_snps[ this_mic - 1 ][ g_num ] ) ) )
        pval_matrix = np.zeros( ( len( valid_snps[ this_mic - 1 ][ g_num ] ), 
                                len( valid_snps[ this_mic - 1 ][ g_num ] ) ) )
        
        for snp1, snp2 in combinations( range( len( valid_snps[ this_mic - 1 ][ g_num ] ) ), 2 ):
            traj1 = valid_snps[ this_mic - 1 ][ g_num ][ snp1 ]
            traj2 = valid_snps[ this_mic - 1 ][ g_num ][ snp2 ]
            pval_matrix[ snp1, snp2 ] = pearsonr( traj1, traj2 )[ 1 ] ** 2
            pval_matrix[ snp2, snp1 ] = pearsonr( traj2, traj1 )[ 1 ] ** 2
            this_matrix[ snp1, snp2 ] = pearsonr( traj1, traj2 )[ 0 ] ** 2
            this_matrix[ snp2, snp1 ] = pearsonr( traj2, traj1 )[ 0 ] ** 2

        for tsnp in range( len( valid_snps[ this_mic - 1 ][ g_num ] ) ):
            ttraj = valid_snps[ this_mic - 1 ][ g_num ][ tsnp ]
            this_matrix[ tsnp, tsnp ] = pearsonr( ttraj, ttraj )[ 0 ] ** 2
            pval_matrix[ tsnp, tsnp ] = pearsonr( ttraj, ttraj )[ 1 ] ** 2
            
        corr_matrix[ this_mic - 1 ][ g_num ] = this_matrix
        sign_matrix[ this_mic - 1 ][ g_num ] = pval_matrix

        # X = this_matrix.copy()
        # X2 = pval_matrix.copy()

        # d = sch.distance.pdist(X)
        # L = sch.linkage(d, method='complete')
        # ind = sch.fcluster(L, 0.5 * d.max(), 'distance')

        # from itertools import product

        # snp_order = np.argsort( ind )
        # list_of_tuples = list( product( snp_order, snp_order ) )

        # corr_X = np.zeros_like( X )
        # corr_p = np.zeros_like( X2 )
        # for idx, (i, j) in enumerate( product( range( len( snp_order ) ), range( len( snp_order ) ) ) ):
        #     corr_X[ i, j ] = X[ list_of_tuples[ idx ][ 0 ], list_of_tuples[ idx ][ 1 ] ]
        #     corr_p[ i, j ] = X2[ list_of_tuples[ idx ][ 0 ], list_of_tuples[ idx ][ 1 ] ]

        # Saving results.
        np.savetxt( '/home/akshitg/lte/metagenomes_lte/freebayes_out/snps/all_alleles/fb_r2_corrs/valid_snps_mic_' 
                    + str( this_mic ) + '_gnm_' + str( g_num ) + '.txt', 
                    valid_snps[ this_mic - 1 ][ g_num ] )
        np.savetxt( '/home/akshitg/lte/metagenomes_lte/freebayes_out/snps/all_alleles/fb_r2_corrs/corr_matrix_mic_' 
                    + str( this_mic ) + '_gnm_' + str( g_num ) + '.txt', 
                    this_matrix )
        np.savetxt( '/home/akshitg/lte/metagenomes_lte/freebayes_out/snps/all_alleles/fb_r2_corrs/sign_matrix_mic_' 
                    + str( this_mic ) + '_gnm_' + str( g_num ) + '.txt', 
                    pval_matrix )

