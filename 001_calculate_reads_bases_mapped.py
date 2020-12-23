import numpy as np
import pandas as pd
import gzip
import gc

col_names_minimap2 = [ str( e ) for e in range( 1, 19 ) ]

def calc_and_write_mapped( mg_num ):
	df = pd.read_csv( 'D19-2600' + str(mg_num) + '_minimap2_default_results.paf', sep='\t', header=None, names=col_names_minimap2 )

	num_total_reads_mapped = len( list( list( df[ '1' ].values ) ) )
	num_unique_reads_mapped = len( set( list( df[ '1' ].values ) ) )

	full = gzip.open( 'D19-2600' + str(mg_num) + '-4151H_combined_mg_reads.fastq.gz', 'r' )
	mg = full.readlines()
	all_reads_in_mg = len( mg ) / 2
	full.close()
	
	if all_reads_in_mg:
		frac_unique_reads_mapped = num_unique_reads_mapped / all_reads_in_mg
		frac_total_reads_mapped = num_total_reads_mapped / all_reads_in_mg
	else:
		frac_unique_reads_mapped = 0.0
		frac_total_reads_mapped = 0.0

	# if df[ '2' ].values:
	frac_bases_mapped_with_gaps = df[ '11' ].values / df[ '2' ].values
	frac_bases_mapped_no_gaps = df[ '10' ].values / df[ '2' ].values
	# else:
	# 	frac_bases_mapped_with_gaps = 0.0
	# 	frac_bases_mapped_no_gaps = 0.0

	map_res_file = open( 'map_res_file.csv', 'a' )
	try:
		res = [ mg_num, num_total_reads_mapped, num_unique_reads_mapped, all_reads_in_mg, frac_total_reads_mapped, frac_unique_reads_mapped , np.mean( frac_bases_mapped_with_gaps ),  np.mean( frac_bases_mapped_no_gaps ) ]
		print(res)
		write_string = ','.join( [ str(e) for e in res ] )
		map_res_file.write( write_string )
		map_res_file.write( '\n' )
	except:
		pass
	map_res_file.close()

 
for mg_num in range(35, 114):
	calc_and_write_mapped( mg_num )
	print( str( mg_num ) + ' DONE.\n')
