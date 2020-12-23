import numpy as np
import pandas as pd

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

for tm00 in allms[-1:]:
    print( 'Processing microcosm ' + str( tm00 ) + '. . .\n' )
    for g_num in g_num_list:
        print( 'Processing genome ' + str( g_num ) + '. . .\n' )
        cumulative_snps = {}
        unique_snps = {}
        is_invalid = False
        for this_idx, this_time in enumerate( tm00 ):
            # Copying over all old SNPs, accumulating them.
            if this_idx > 0:
                cumulative_snps[ this_time ] = cumulative_snps[ tm00[ this_idx - 1 ] ][ : ]
            else:
                cumulative_snps[ this_time ] = []

            # Going from time point to time point.
            try:
                df = pd.read_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/snps/fb_processed/mg_' 
                                  + str( this_time ) + '_gnm_' + str( g_num ) + 
                                  '_processed_snps.tsv', sep='\t', index_col=0 )

                if len( df ) == 0:
                    continue
            except:
                continue

            unique_snps[ this_time ] = []
            # Now go SNP by SNP.
            for idx, row in  df.iterrows():
                tsnp = str( row[ 'contig' ] ) + '_' + str( row[ 'gene_num' ] ) + '_' + str( row[ 'pos' ] ) + '_' + str( row[ 'ref_base' ] ) + '_' + str( row[ 'alt_base' ] )
                if this_idx == 0:
                    cumulative_snps[ this_time ].append( tsnp )
                    unique_snps[ this_time ].append( tsnp )
                elif tsnp not in cumulative_snps[ tm00[ this_idx - 1 ] ]:
                    cumulative_snps[ this_time ].append( tsnp )
                    unique_snps[ this_time ].append( tsnp )

        for this_idx, this_time in enumerate( tm00 ):
            try:
                # Going from time point to time point.
                df = pd.read_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/snps/fb_processed/mg_' 
                                  + str( this_time ) + '_gnm_' + str( g_num ) + 
                                  '_processed_snps.tsv', sep='\t', index_col=0 )

                if len( df ) == 0:
                    is_invalid = True
                    continue
            except:
                is_invalid = True
                continue

            ndf = pd.DataFrame()

            # Now go SNP by SNP.
            for idx, row in  df.iterrows():
                if not len( df ):
                    break
                tsnp = str( row[ 'contig' ] ) + '_' + str( row[ 'gene_num' ] ) + '_' + str( row[ 'pos' ] ) + '_' + str( row[ 'ref_base' ] ) + '_' + str( row[ 'alt_base' ] )
                if tsnp in unique_snps[ this_time ]:
                    ndf = ndf.append( df.loc[ idx ] )
            if len( ndf ):
                ndf = ndf[ df.columns.tolist() ]
                ndf[ 'gene_num' ] = ndf[ 'gene_num' ].astype( int )
                ndf[ 'pos' ] = ndf[ 'pos' ].astype( int )
                ndf.to_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/snps/fb_unique_no_threshold/mg_' 
                            + str( this_time ) + '_gnm_' + str( g_num ) + '_unique_snps.tsv', sep='\t' )

