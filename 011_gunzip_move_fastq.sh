for d in D19-*/ ; do
    for f in ~/lte/metagenomes_lte/"${d}"/*.gz; do
	STEM=$(basename "${f}" .gz)
	mkdir ~/corderolab/akshitg/"${d}"
	gunzip -c "${f}" > ~/corderolab/akshitg/"${d}"/"${STEM}"
    echo "${d} done."
    done
done

#for f in ~/lte/metagenomes_lte/D19-260037-4151H/*.gz; do
#   STEM=$(basename "${f}" .gz)
#   gunzip -c "${f}" > ~/corderolab/akshitg/D19-260037-4151H/"${STEM}"
#   echo "DONE ${f}"
#done
