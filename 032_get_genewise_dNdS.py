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

dNdS_gene_df = pd.DataFrame()

microcosm_list, metagenome_list, genome_list, gene_list, gene_code_list, ko_list = [], [], [], [], [], []
module_list, map_list, snp_num_list = [], [], []
S_list, sd_list, N_list, nd_list, pNpS_list, dNdS_list = [], [], [], [], [], []

this_time = int( sys.argv[ 1 ] )

for this_mic in mic_list:
    if this_time in eval( this_mic ):
        break

SNP_THRESHOLD = 5

for g_num in g_num_list:
    df = pd.read_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/snps/genewise_snps/mg_' 
                      + str( this_time ) + 
                      '_snp_gene_table.tsv', sep='\t', index_col=0 ).fillna( '' )

    genome_df = df[ df[ 'genome' ] == g_num ]

    # Getting a list of mutated genes for this genome.
    this_mutated_genes = sorted( list( set( genome_df[ 'gene_num' ].values ) ) )

    for this_gene in this_mutated_genes:
        S, N, nd, sd = 0.0, 0.0, 0.0, 0.0
        these_gene_rows = genome_df[ genome_df[ 'gene_num' ] == this_gene ]
        if len( these_gene_rows ) > SNP_THRESHOLD:
            # Now go SNP by SNP.
            for idx, row in  these_gene_rows.iterrows():
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

            # Now calculating the gene-wide pN/pS and dN/dS ratios for this gene.
            if sd:
                this_pNpS = ( nd / N ) / ( sd / S )
                if ( 1 - ( 4 / 3 ) * ( sd / S ) ) <= 0.0 or ( 1 - ( 4 / 3 ) * ( nd / N ) ) <= 0.0:
                    this_dNdS = 0.0
                else:
                    this_dNdS = np.log( 1 - ( 4 / 3 ) * ( nd / N ) ) / np.log( 1 - ( 4 / 3 ) * ( sd / S ) )
            else:
                this_pNpS = 0.0
                this_dNdS = 0.0

            microcosm_list.append( this_mic )
            metagenome_list.append( this_time )
            genome_list.append( g_num )
            gene_list.append( this_gene )
            gene_code_list.append( row[ 'gene_code' ] )
            ko_list.append( row[ 'ko' ] )
            module_list.append( row[ 'module' ] )
            map_list.append( row[ 'map' ] )
            snp_num_list.append( len( these_gene_rows ) )
            S_list.append( S )
            sd_list.append( sd )
            N_list.append( N )
            nd_list.append( nd )
            pNpS_list.append( this_pNpS )
            dNdS_list.append( this_dNdS )

dNdS_gene_df[ 'microcosm' ] = microcosm_list
dNdS_gene_df[ 'metagenome' ] = metagenome_list
dNdS_gene_df[ 'genome' ] = genome_list
dNdS_gene_df[ 'gene_num' ] = gene_list
dNdS_gene_df[ 'gene_code' ] = gene_code_list
dNdS_gene_df[ 'ko' ] = ko_list
dNdS_gene_df[ 'module' ] = module_list
dNdS_gene_df[ 'map' ] = map_list
dNdS_gene_df[ 'num_snps' ] = snp_num_list
dNdS_gene_df[ 'S' ] = S_list
dNdS_gene_df[ 'sd' ] = sd_list
dNdS_gene_df[ 'N' ] = N_list
dNdS_gene_df[ 'nd' ] = nd_list
dNdS_gene_df[ 'pNpS' ] = pNpS_list
dNdS_gene_df[ 'dNdS' ] = dNdS_list

dNdS_gene_df.to_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/mg_' + 
                    str( this_time ) + '_dNdS_gene_table.tsv', sep='\t' )

