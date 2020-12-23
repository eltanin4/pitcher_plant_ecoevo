import numpy as np
import pandas as pd

###############################################################################
# Saving all synonymous reference sites.
###############################################################################
syn_sites = { 
    'ATA': 0.67, 'ATC': 0.67, 'ATT': 0.67, 'ATG': 0.00, 
    'ACA': 1.00, 'ACC': 1.00, 'ACG': 1.00, 'ACT': 1.00, 
    'AAC': 0.33, 'AAT': 0.33, 'AAA': 0.33, 'AAG': 0.33, 
    'AGC': 0.33, 'AGT': 0.33, 'AGA': 0.83, 'AGG': 0.67, 
    'CTA': 1.33, 'CTC': 1.00, 'CTG': 1.33, 'CTT': 1.00, 
    'CCA': 1.00, 'CCC': 1.00, 'CCG': 1.00, 'CCT': 1.00, 
    'CAC': 0.33, 'CAT': 0.33, 'CAA': 0.33, 'CAG': 0.33, 
    'CGA': 1.50, 'CGC': 1.00, 'CGG': 1.33, 'CGT': 1.00, 
    'GTA': 1.00, 'GTC': 1.00, 'GTG': 1.00, 'GTT': 1.00, 
    'GCA': 1.00, 'GCC': 1.00, 'GCG': 1.00, 'GCT': 1.00, 
    'GAC': 0.33, 'GAT': 0.33, 'GAA': 0.33, 'GAG': 0.33, 
    'GGA': 1.00, 'GGC': 1.00, 'GGG': 1.00, 'GGT': 1.00, 
    'TCA': 1.00, 'TCC': 1.00, 'TCG': 1.00, 'TCT': 1.00, 
    'TTC': 0.33, 'TTT': 0.33, 'TTA': 0.67, 'TTG': 0.67, 
    'TAC': 1.00, 'TAT': 1.00, 'TAA': 0.00, 'TAG': 0.00, 
    'TGC': 0.50, 'TGT': 0.50, 'TGA': 0.00, 'TGG': 0.00,
} 

###############################################################################
# Building a library of microcosm to metagenome numbers.
###############################################################################
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

###############################################################################
# Building a dictionary to get microcosm number.
###############################################################################
mic_dict = { '[35, 43, 51, 59, 67, 75, 83, 91]' : 'M01',
             '[36, 44, 52, 60, 68, 76, 84, 92]' : 'M02',
             '[37, 45, 53, 61, 69, 77, 85, 93]' : 'M03',
             '[38, 46, 54, 62, 70, 78, 86, 94]' : 'M04',
             '[39, 47, 55, 63, 71, 79, 87, 95]' : 'M05',
             '[40, 48, 56, 64, 72, 80, 88, 96]' : 'M06',
             '[41, 49, 57, 65, 73, 81, 89, 97]' : 'M07',
             '[42, 50, 58, 66, 74, 82, 90, 98]' : 'M08',
             '[99, 101, 103, 105, 107, 109, 111, 113]' : 'M09',
             '[100, 102, 104, 106, 108, 110, 112, 114]' : 'M10'
           }

###############################################################################
# Setting up key variables, and initialising them.
###############################################################################
allms = [ m01, m02, m03, m04, m05, m06, m07, m08, m09, m10 ]
g_num_list = list( range( 1, 25 ) ) + list( range( 26, 35 ) )

dNdS_df = pd.DataFrame()
microcosm_list, metagenome_list, genome_list, S_list, sd_list, N_list, nd_list, pNpS_list, dNdS_list = [], [], [], [], [], [], [], [], []

###############################################################################
# Calculating pN/pS and dN/dS genome by genome.
###############################################################################
for tm00 in allms:
    print( 'Processing microcosm ' + str( tm00 ) + '. . .\n' )
    for g_num in g_num_list:
        print( 'Processing genome ' + str( g_num ) + '. . .\n' )
        for this_idx, this_time in enumerate( tm00 ):
            # Going from time point to time point.
            try:
                df = pd.read_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/snps/fb_unique/mg_' 
                                  + str( this_time ) + '_gnm_' + str( g_num ) + 
                                  '_unique_snps.tsv', sep='\t', index_col=0 )

                if len( df ) == 0:
                    continue
            except:
                continue

            S, N, nd, sd = 0.0, 0.0, 0.0, 0.0
            # Now go SNP by SNP.
            for idx, row in  df.iterrows():
                # Excluding stop codon mutations from the analysis.
                if row[ 'snp_type' ] == 'stop':
                    continue

                # Collecting SNP and reference codon and base to calculate the number of expected sites.
                this_ref_codon = row[ 'ref_codon' ]
                this_ref_base = row[ 'ref_base' ]

                S += syn_sites[ this_ref_codon ]
                N += 3.0 - syn_sites[ this_ref_codon ]

                if row[ 'snp_type' ] == 'syn':
                    sd += 1.0
                elif row[ 'snp_type' ] == 'non':
                    nd += 1.0

            # Now calculating the genome-wide pN/pS and dN/dS ratios for this genome.
            if sd:
                this_pNpS = ( nd / N ) / ( sd / S )
                this_dNdS = np.log( 1 - ( 4 / 3 ) * ( nd / N ) ) / np.log( 1 - ( 4 / 3 ) * ( sd / S ) )
            else:
                this_pNpS = 0.0
                this_dNdS = 0.0
            microcosm_list.append( mic_dict[ str( tm00 ) ] )
            metagenome_list.append( this_time )
            genome_list.append( g_num )
            S_list.append( S )
            sd_list.append( sd )
            N_list.append( N )
            nd_list.append( nd )
            pNpS_list.append( this_pNpS )
            dNdS_list.append( this_dNdS )

###############################################################################
# Collecting results into a dataframe.
###############################################################################
dNdS_df[ 'microcosm' ] = microcosm_list
dNdS_df[ 'metagenome' ] = metagenome_list
dNdS_df[ 'genome_num' ] = genome_list
dNdS_df[ 'S' ] = S_list
dNdS_df[ 'sd' ] = sd_list
dNdS_df[ 'N' ] = N_list
dNdS_df[ 'nd' ] = nd_list
dNdS_df[ 'pNpS' ] = pNpS_list
dNdS_df[ 'dNdS' ] = dNdS_list

###############################################################################
# Saving the processed dataframe with results.
###############################################################################
dNdS_df.to_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/fb_dnds_table.tsv', sep='\t' )
print( 'Final processed table saved.' )

