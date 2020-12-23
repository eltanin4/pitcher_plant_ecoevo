import numpy as np
import pandas as pd
import gzip
import sys

mg_num = int( sys.argv[ 1 ] )
g_num = int( sys.argv[ 2 ] )

df = pd.read_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/bt2_snps/tsvs/mg_' 
                  + str( mg_num ) + '_gnm_' + str( g_num ) + '_only_snps_par.tsv', sep='\t' )

PVAL_THRESH = 30
tgenedf = df.copy()
filtdf = tgenedf[ np.logical_and( tgenedf[ 'QUAL' ] >= PVAL_THRESH, tgenedf[ 'DP' ] >= 10 ) ]
filtdf = filtdf[ filtdf[ 'TYPE' ] == 'snp' ]

###############################################################################
# Reading genes
###############################################################################
genes_file = open( '/home/akshitg/main/lte/sags_lte/D19-2600' + ( 2 - len( str( g_num ) ) ) * '0' +  str( g_num ) + '/genes.faa', 'r' )
gene_names_raw = genes_file.readlines()
gene_content_raw = ''.join( gene_names_raw ).split( '>' )[ 1: ]

gene_dict = {}
for gene_idx in range( len( gene_content_raw ) ):
    gene_id = gene_content_raw[ gene_idx ].split( '\n' )[ 0 ]
    gene_seq = ''.join( gene_content_raw[ gene_idx ].split( '\n' )[ 1: ] )
    gene_dict[ gene_id ] = gene_seq

gene_df = pd.DataFrame()
gene_nums, gene_seqs, gene_starts, gene_ends, gene_contigs = [], [], [], [], []
gene_names = {}
for gene_idx, gene_name in enumerate( gene_dict ):
    gene_nums.append( gene_idx + 1 )
    gene_names[ gene_idx + 1 ] = gene_name
    gene_seqs.append( gene_dict[ gene_name ] )
    gene_starts.append( int( gene_name.split( ' # ' )[1] ) )
    gene_ends.append( int( gene_name.split( ' # ' )[2] ) )
    gene_contigs.append( '_'.join( gene_name.split( ' # ' )[0].split( '_' )[:-1] ) )
gene_df[ 'gene_num' ] = gene_nums
gene_df[ 'seq' ] = gene_seqs
gene_df[ 'start' ] = gene_starts
gene_df[ 'end' ] = gene_ends
gene_df[ 'contig' ] = gene_contigs

###############################################################################
# Reading proteins
###############################################################################
proteins_file = open( '/home/akshitg/main/lte/sags_lte/D19-2600' + ( 2 - len( str( g_num ) ) ) * '0' +  str( g_num ) + '/proteins.faa', 'r' )
protein_names_raw = proteins_file.readlines()
protein_content_raw = ''.join( protein_names_raw ).split( '>' )[ 1: ]

protein_dict = {}
for gene_idx in range( len( protein_content_raw ) ):
    gene_id = gene_idx + 1
    aa_seq = ''.join( protein_content_raw[ gene_idx ].split( '\n' )[ 1: ] )
    protein_dict[ gene_id ] = aa_seq

###############################################################################
# Getting the protein table
###############################################################################
from collections import Counter, defaultdict

table = { 
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
} 

###############################################################################
# Collecting all information and assigning SNPs
###############################################################################
curr_gene_num = 1
curr_start = gene_df[ gene_df[ 'gene_num' ] == curr_gene_num ].iloc[ 0, 2 ]
curr_end = gene_df[ gene_df[ 'gene_num' ] == curr_gene_num ].iloc[ 0, 3 ]

snp_df = pd.DataFrame()
snp_contigs, snp_gnums, snp_poses, snp_refs, snp_alts, snp_aarefs, snp_aaalts, snp_types = [], [], [], [], [], [], [], []
snp_refcounts, snp_altcounts, snp_depths, snp_reffreqs, snp_altfreqs, snp_pis = [], [], [], [], [], []
snp_ref_cods, snp_alt_cods = [], []

for idx, row in  filtdf.iterrows():
    not_genic_snp = False
    rel_gene_df = gene_df[ gene_df[ 'contig' ] == row[ '#CHROM' ] ]
    
    if not len( rel_gene_df ):
        continue
    
    curr_gene_num = rel_gene_df.iloc[ 0, 0 ]
    curr_start = rel_gene_df[ rel_gene_df[ 'gene_num' ] == curr_gene_num ].iloc[ 0, 2 ]
    curr_end = rel_gene_df[ rel_gene_df[ 'gene_num' ] == curr_gene_num ].iloc[ 0, 3 ]
    
    if row[ 'POS' ] < curr_start:
        not_genic_snp = True
    
    
    if gene_dict[ gene_names[ curr_gene_num ] ][ int( ( row[ 'POS' ] - curr_start ) / 3 ) * 3 : int( ( row[ 'POS' ] - curr_start ) / 3 ) * 3 + 3 ] == 'NGG':
       not_genic_snp = True
    
    while row[ 'POS' ] > curr_end and not not_genic_snp:
        if curr_end == rel_gene_df.iloc[ -1, 3]:
            not_genic_snp = True
            break
        curr_gene_num += 1
        curr_start = rel_gene_df[ rel_gene_df[ 'gene_num' ] == curr_gene_num ].iloc[ 0, 2 ]
        curr_end = rel_gene_df[ rel_gene_df[ 'gene_num' ] == curr_gene_num ].iloc[ 0, 3 ]
        if curr_start > row[ 'POS' ]:
            not_genic_snp = True
            break
    
    
    if len( row[ 'REF' ] ) == 1 and not not_genic_snp:
        snp_contigs.append( row[ '#CHROM' ] )
        snp_gnums.append( rel_gene_df[ rel_gene_df[ 'gene_num' ] == curr_gene_num ].iloc[ 0, 0 ] )
        snp_poses.append( row[ 'POS' ] )
        snp_refs.append( row[ 'REF' ] )
        snp_alts.append( row[ 'ALT' ] )
        try:
            snp_aarefs.append( protein_dict[ curr_gene_num ][ int( ( row[ 'POS' ] - curr_start ) / 3 ) ] )
        except:
            print( curr_gene_num, row[ 'POS' ], curr_start, curr_end )

        # Reference codon
        ref_cod = gene_dict[ gene_names[ curr_gene_num ] ][ int( ( row[ 'POS' ] - curr_start ) / 3 ) * 3 : int( ( row[ 'POS' ] - curr_start ) / 3 ) * 3 + 3 ]
        if not table[ ref_cod ] == protein_dict[ curr_gene_num ][ int( ( row[ 'POS' ] - curr_start ) / 3 ) ]:
            if not protein_dict[ curr_gene_num ][ int( ( row[ 'POS' ] - curr_start ) / 3 ) ] == 'M':            
                print( table[ ref_cod ] )
                print( protein_dict[ curr_gene_num ][ int( ( row[ 'POS' ] - curr_start ) / 3 ) ] )
                print( curr_gene_num,
                       row[ '#CHROM' ],
                       row[ 'POS' ] )
                raise ValueError( 'Mismatch between expected codon and observed codon at reference' )
        
        # Altered codon
        gene_dict[ gene_names[ curr_gene_num ] ][ int( ( row[ 'POS' ] - curr_start ) / 3 ) * 3 : int( ( row[ 'POS' ] - curr_start ) / 3 ) * 3 + 3 ]
        alt_gene_seq = gene_dict[ gene_names[ curr_gene_num ] ][:]
        alt_gene_seq = alt_gene_seq[ : row[ 'POS' ] - curr_start ] + row[ 'ALT' ]  + alt_gene_seq[ row[ 'POS' ] - curr_start + 1 : ]
        alt_cod = alt_gene_seq[ int( ( row[ 'POS' ] - curr_start ) / 3 ) * 3 : int( ( row[ 'POS' ] - curr_start ) / 3 ) * 3 + 3 ]
        snp_aaalts.append( table[ alt_cod ] )

        # Writing the reference and alternate codons.
        snp_ref_cods.append( ref_cod )
        snp_alt_cods.append( alt_cod )

        # Calculating allele frequencies and at-site nucleotide diversity.
        p_ref = float( int( row[ 'RO' ] ) / int( row[ 'DP' ] ) )
        p_alt = float( int( row[ 'AO' ] ) / int( row[ 'DP' ] ) )
        pi_site = 1 - ( p_ref ** 2 ) - ( p_alt ** 2 )

        # Writing the counts, depths and frequencies.
        snp_refcounts.append( int( row[ 'RO' ] ) )
        snp_altcounts.append( int( row[ 'AO' ] ) )
        snp_depths.append( int( row[ 'DP' ] ) )
        snp_reffreqs.append( p_ref )
        snp_altfreqs.append( p_alt )
        snp_pis.append( pi_site )


        if protein_dict[ curr_gene_num ][ int( ( row[ 'POS' ] - curr_start ) / 3 ) ] == table[ alt_cod ]:
            snp_types.append( 'syn' )
        elif table[ alt_cod ] == '*':
            snp_types.append( 'stop' )
        else:
            snp_types.append( 'non' )

print( 'Finished processing SNPs for genome ' + str( g_num ) + '.' )

snp_df[ 'contig' ] = snp_contigs
snp_df[ 'gene_num' ] = snp_gnums
snp_df[ 'pos' ] = snp_poses
snp_df[ 'ref_base' ] = snp_refs
snp_df[ 'alt_base' ] = snp_alts
snp_df[ 'ref_codon' ] = snp_ref_cods
snp_df[ 'alt_codon' ] = snp_alt_cods
snp_df[ 'ref_aa' ] = snp_aarefs
snp_df[ 'alt_aa' ] = snp_aaalts
snp_df[ 'ref_count' ] = snp_refcounts
snp_df[ 'alt_count' ] = snp_altcounts
snp_df[ 'locus_depth' ] = snp_depths
snp_df[ 'ref_freq' ] = snp_reffreqs
snp_df[ 'alt_freq' ] = snp_altfreqs
snp_df[ 'pi_site' ] = snp_pis
snp_df[ 'snp_type' ] = snp_types

snp_df.to_csv( '/home/akshitg/lte/metagenomes_lte/freebayes_out/snps/bt2/all_alleles/fb_processed/mg_' 
               + str( mg_num ) + '_gnm_' + str( g_num ) + '_processed_snps.tsv', sep='\t' )
