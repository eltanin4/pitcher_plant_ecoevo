import numpy as np
import pandas as pd
import sys

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

tm00 = allms[ int( sys.argv[ 1 ] ) - 1 ]
print( 'Processing microcosm ' + str( tm00 ) )

for g_num in g_num_list:
    print( 'Processing genome ' + str( g_num ) + '.' )
    
    snp_detects = {}
    snp_detectpoints = {}

    is_invalid = False
    for this_idx, this_time in enumerate( tm00 ):
        # Going from time point to time point.
        try:
            df = pd.read_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/snps/all_alleles/fb_processed/mg_'
                              + str( this_time ) + '_gnm_' + str( g_num ) +
                              '_processed_snps.tsv', sep='\t', index_col=0 )

            if len( df ) == 0:
                is_invalid = True
                continue
        except:
            is_invalid = True
            continue

        # Now go SNP by SNP.
        for idx, row in  df.iterrows():
            tsnp = str( row[ 'contig' ] ) + '_' + str( row[ 'gene_num' ] ) + '_' + str( row[ 'pos' ] ) + '_' + str( row[ 'ref_base' ] ) + '_' + str( row[ 'alt_base' ] )
            if tsnp in snp_detects:
                snp_detects[ tsnp ] += 1
            else:
                snp_detects[ tsnp ] = 1
                snp_detectpoints[ tsnp ] = this_idx

    if is_invalid:
        continue

    # Choosing SNP threshold stringently. I also want long-lived trajectories if I can get them.
    # Importantly, I want to purge trajectories that die quickly.
    THRESHOLD_N = 2
    
    # Reducing the treshold for microcosms that are missing one time point, M06 (40) and M10 (106).
    if int( sys.argv[ 1 ] ) in [ 6, 10 ]:
        THRESHOLD_N -= 1
    bad_snps = []
    for tsnp in snp_detects:
        if snp_detectpoints[ tsnp ] >= len( tm00 ) - THRESHOLD_N:        
            if snp_detects[ tsnp ] < len( tm00 ) - snp_detectpoints[ tsnp ]:
                bad_snps.append( tsnp )
        else:
            if snp_detects[ tsnp ] < THRESHOLD_N: 
                bad_snps.append( tsnp )

    for tsnp in bad_snps:
        snp_detects.pop( tsnp, None )
        snp_detectpoints.pop( tsnp, None )

    tdf = {}
    for this_idx, this_time in enumerate( tm00 ):
        try:
            # Going from time point to time point.
            tdf[ this_idx ] = pd.read_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/snps/all_alleles/fb_processed/mg_'
                                           + str( this_time ) + '_gnm_' + str( g_num ) +
                                           '_processed_snps.tsv', sep='\t', index_col=0 )

            if len( df ) == 0:
                is_invalid = True
                continue
        except:
            is_invalid = True
            continue

    snp_traj_df = pd.DataFrame()
    snp_contigs, snp_gnums, snp_poses, snp_refs, snp_alts = [], [], [], [], []
    snp_ref_cods, snp_alt_cods, snp_aarefs, snp_aaalts, snp_types = [], [], [], [], []

    t1_refs, t2_refs, t3_refs, t4_refs, t5_refs, t6_refs, t7_refs, t8_refs = [], [], [], [], [], [], [], []
    t1_alts, t2_alts, t3_alts, t4_alts, t5_alts, t6_alts, t7_alts, t8_alts = [], [], [], [], [], [], [], []
    t1_depths, t2_depths, t3_depths, t4_depths, t5_depths, t6_depths, t7_depths, t8_depths = [], [], [], [], [], [], [], []

    # Now go SNP by SNP.
    for tsnp in snp_detects:

        rel_df = tdf[ snp_detectpoints[ tsnp ] ]
        rel_row = rel_df[ np.logical_and.reduce( [ ( rel_df[ 'contig' ] == '_'.join( tsnp.split( '_' )[ :-4 ] ) ).values, 
                                                   ( rel_df[ 'gene_num' ] == int( tsnp.split( '_' )[ -4 ] ) ).values,
                                                   ( rel_df[ 'pos' ] == int( tsnp.split( '_' )[ -3 ] ) ).values,
                                                   ( rel_df[ 'ref_base' ] == tsnp.split( '_' )[ -2 ] ).values,
                                                   ( rel_df[ 'alt_base' ] == tsnp.split( '_' )[ -1 ] ).values ] ) ]
        
        snp_contigs.append( rel_row[ 'contig' ].values[0] )
        snp_gnums.append( int( rel_row[ 'gene_num' ].values[0] ) )
        snp_poses.append( int( rel_row[ 'pos' ].values[0] ) )
        snp_refs.append( rel_row[ 'ref_base' ].values[0] )
        snp_alts.append( rel_row[ 'alt_base' ].values[0] )
        snp_ref_cods.append( rel_row[ 'ref_codon' ].values[0] )
        snp_alt_cods.append( rel_row[ 'alt_codon' ].values[0] )
        snp_aarefs.append( rel_row[ 'ref_aa' ].values[0] )
        snp_aaalts.append( rel_row[ 'alt_aa' ].values[0] )
        snp_types.append( rel_row[ 'snp_type' ].values[0] )

        for this_idx, this_time in enumerate( tm00 ):
            if this_idx < snp_detectpoints[ tsnp ]:
                eval( 't' + str( this_idx + 1 ) + '_refs.append( 0.0 )' )
                eval( 't' + str( this_idx + 1 ) + '_alts.append( 0.0 )' )
                eval( 't' + str( this_idx + 1 ) + '_depths.append( 0 )' )
            elif this_idx == snp_detectpoints[ tsnp ]:
                eval( 't' + str( this_idx + 1 ) + '_refs.append( float( rel_row[ \'ref_freq\' ].values[0] ) )' )
                eval( 't' + str( this_idx + 1 ) + '_alts.append( float( rel_row[ \'alt_freq\' ].values[0] ) )' )
                eval( 't' + str( this_idx + 1 ) + '_depths.append( int( rel_row[ \'locus_depth\' ].values[0] ) )' )
            else:
                rel_df = tdf[ this_idx ]
                rel_row = rel_df[ np.logical_and.reduce( [ ( rel_df[ 'contig' ] == '_'.join( tsnp.split( '_' )[ :-4 ] ) ).values, 
                                           ( rel_df[ 'gene_num' ] == int( tsnp.split( '_' )[ -4 ] ) ).values,
                                           ( rel_df[ 'pos' ] == int( tsnp.split( '_' )[ -3 ] ) ).values,
                                           ( rel_df[ 'ref_base' ] == tsnp.split( '_' )[ -2 ] ).values,
                                           ( rel_df[ 'alt_base' ] == tsnp.split( '_' )[ -1 ] ).values ] ) ]
                if len( rel_row ) == 0:
                    eval( 't' + str( this_idx + 1 ) + '_refs.append( 0.0 )' )
                    eval( 't' + str( this_idx + 1 ) + '_alts.append( 0.0 )' )
                    eval( 't' + str( this_idx + 1 ) + '_depths.append( 0 )' )
                else:
                    eval( 't' + str( this_idx + 1 ) + '_refs.append( float( rel_row[ \'ref_freq\' ].values[0] ) )' )
                    eval( 't' + str( this_idx + 1 ) + '_alts.append( float( rel_row[ \'alt_freq\' ].values[0] ) )' )
                    eval( 't' + str( this_idx + 1 ) + '_depths.append( int( rel_row[ \'locus_depth\' ].values[0] ) )' )

    snp_traj_df[ 'contig' ] = snp_contigs
    snp_traj_df[ 'gene_num' ] = snp_gnums
    snp_traj_df[ 'pos' ] = snp_poses
    snp_traj_df[ 'ref_base' ] = snp_refs
    snp_traj_df[ 'alt_base' ] = snp_alts
    snp_traj_df[ 'ref_codon' ] = snp_ref_cods
    snp_traj_df[ 'alt_codon' ] = snp_alt_cods
    snp_traj_df[ 'ref_aa' ] = snp_aarefs
    snp_traj_df[ 'alt_aa' ] = snp_aaalts
    snp_traj_df[ 'snp_type' ] = snp_types

    snp_traj_df[ 't1_ref_freq' ] = t1_refs
    snp_traj_df[ 't2_ref_freq' ] = t2_refs
    snp_traj_df[ 't3_ref_freq' ] = t3_refs
    snp_traj_df[ 't4_ref_freq' ] = t4_refs
    snp_traj_df[ 't5_ref_freq' ] = t5_refs
    snp_traj_df[ 't6_ref_freq' ] = t6_refs
    snp_traj_df[ 't7_ref_freq' ] = t7_refs
    snp_traj_df[ 't8_ref_freq' ] = t8_refs

    snp_traj_df[ 't1_alt_freq' ] = t1_alts
    snp_traj_df[ 't2_alt_freq' ] = t2_alts
    snp_traj_df[ 't3_alt_freq' ] = t3_alts
    snp_traj_df[ 't4_alt_freq' ] = t4_alts
    snp_traj_df[ 't5_alt_freq' ] = t5_alts
    snp_traj_df[ 't6_alt_freq' ] = t6_alts
    snp_traj_df[ 't7_alt_freq' ] = t7_alts
    snp_traj_df[ 't8_alt_freq' ] = t8_alts

    snp_traj_df[ 't1_locus_depth' ] = t1_depths
    snp_traj_df[ 't2_locus_depth' ] = t2_depths
    snp_traj_df[ 't3_locus_depth' ] = t3_depths
    snp_traj_df[ 't4_locus_depth' ] = t4_depths
    snp_traj_df[ 't5_locus_depth' ] = t5_depths
    snp_traj_df[ 't6_locus_depth' ] = t6_depths
    snp_traj_df[ 't7_locus_depth' ] = t7_depths
    snp_traj_df[ 't8_locus_depth' ] = t8_depths

    snp_traj_df.to_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/snps/all_alleles/fb_trajs/mic_'
                        + sys.argv[ 1 ] + '_gnm_' + str( g_num ) + '_snp_trajectories.tsv', sep='\t' )

