import numpy as np
import pandas as pd
import sys

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
M01 = [ 35, 43, 51, 59, 67, 75, 83, 91 ]
M02 = [ 36, 44, 52, 60, 68, 76, 84, 92 ]
M03 = [ 37, 45, 53, 61, 69, 77, 85, 93 ]
M04 = [ 38, 46, 54, 62, 70, 78, 86, 94 ]
M05 = [ 39, 47, 55, 63, 71, 79, 87, 95 ]
M06 = [ 40, 48, 56, 64, 72, 80, 88, 96 ]
M07 = [ 41, 49, 57, 65, 73, 81, 89, 97 ]
M08 = [ 42, 50, 58, 66, 74, 82, 90, 98 ]
M09 = [ 99, 101, 103, 105, 107, 109, 111, 113 ]
M10 = [ 100, 102, 104, 106, 108, 110, 112, 114 ]

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
kodf = pd.read_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/gene_to_ko_ref.tsv', 
                    sep='\t', index_col=0 ).fillna('')

g_num_list = list( range( 1, 25 ) ) + list( range( 26, 35 ) )

mic_list = [ 'M01', 'M02', 'M03', 'M04', 'M05', 'M06', 'M07', 'M08', 'M09', 'M10' ]

snp_gene_df = pd.DataFrame()

microcosm_list, metagenome_list, genome_list, gene_list, gene_code_list, ko_list = [], [], [], [], [], []
module_list, map_list, snp_list = [], [], []
snp_refs, snp_alts, snp_ref_cods, snp_alt_cods, snp_aarefs, snp_aaalts = [], [], [], [], [], []

this_time = int( sys.argv[ 1 ] )

for this_mic in mic_list:
    if this_time in eval( this_mic ):
        break

for g_num in g_num_list:
    try:
        df = pd.read_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/snps/fb_unique_no_threshold/mg_' 
                          + str( this_time ) + '_gnm_' + str( g_num ) + 
                          '_unique_snps.tsv', sep='\t', index_col=0 )

        if len( df ) == 0:
            continue
    except:
        continue

    # akshit_gene_code represents the gene numbers I mapped from the Prodigal detections.
    # Each genome has a unique gene number, the last number being the number of detected genes.
    # A subset of these are annotated by eggNOG, whose codes are in eggnog_gene_code.
    for idx, row in  df.iterrows():
        this_gene_code = str( df.loc[ idx ][ 'contig' ] ) + '_' + str( df.loc[ idx ][ 'gene_num' ] )
        if kodf[ kodf[ 'akshit_gene_code' ] == this_gene_code ][ 'ko' ].values:
            microcosm_list.append( this_mic )
            metagenome_list.append( this_time )
            genome_list.append( g_num )
            gene_list.append( df.loc[ idx ][ 'gene_num' ] )
            gene_code_list.append( this_gene_code )
            ko_list.append( kodf[ kodf[ 'akshit_gene_code' ] == this_gene_code ][ 'ko' ].values[0] )
            module_list.append( kodf[ kodf[ 'akshit_gene_code' ] == this_gene_code ][ 'module' ].values[0] )
            map_list.append( kodf[ kodf[ 'akshit_gene_code' ] == this_gene_code ][ 'map' ].values[0] )
            snp_list.append( row[ 'snp_type' ] )
            snp_refs.append( row[ 'ref_base' ] )
            snp_alts.append( row[ 'alt_base' ] )
            snp_ref_cods.append( row[ 'ref_codon' ] )
            snp_alt_cods.append( row[ 'alt_codon' ] )
            snp_aarefs.append( row[ 'ref_aa' ] )
            snp_aaalts.append( row[ 'alt_aa' ] )

snp_gene_df[ 'microcosm' ] = microcosm_list
snp_gene_df[ 'metagenome' ] = metagenome_list
snp_gene_df[ 'genome' ] = genome_list
snp_gene_df[ 'gene_num' ] = gene_list
snp_gene_df[ 'gene_code' ] = gene_code_list
snp_gene_df[ 'ko' ] = ko_list
snp_gene_df[ 'module' ] = module_list
snp_gene_df[ 'map' ] = map_list
snp_gene_df[ 'ref_base' ] = snp_refs
snp_gene_df[ 'alt_base' ] = snp_alts
snp_gene_df[ 'ref_codon' ] = snp_ref_cods
snp_gene_df[ 'alt_codon' ] = snp_alt_cods
snp_gene_df[ 'ref_aa' ] = snp_aarefs
snp_gene_df[ 'alt_aa' ] = snp_aaalts
snp_gene_df[ 'snp_type' ] = snp_list

snp_gene_df.to_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/mg_' + 
                    str( this_time ) + '_snp_gene_table.tsv', sep='\t' )

