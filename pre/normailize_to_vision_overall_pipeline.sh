###### read parameters from inputs
script_dir=$1
working_dir=$2
input_dir=$3
input_file_list=$4
overall_upper=$5
overall_lower=$6
select_method=$7
user_given_global_ref=$8
user_given_mark_ref_list=$9

script_dir='/storage/home/g/gzx103/group/projects/vision/normalize_to_vision/src/'
bw_folder='/storage/home/g/gzx103/group/projects/vision/normalize_to_vision/signals/'
bedfile_200bp='mm10.200bin.rand.sort.bed'
bw_file_list='bw_file_list.txt'
bw_file_allref='ideasVisionV20p8NormImkK27ac.bw'
upperlim1=100
lowerlim1=0
overall_ref=''

### get random windows
bedtools makewindows -g ~/group/genome/mm10/mm10.chrom.sizes -w 200 > mm10.200bin.bed
Rscript get_random_200bp_bins.R
cat mm10.200bin.rand.bed | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' | sort -k1,1 -k2,2n > $bedfile_200bp

### save output bed 
sort -k4,4 $bedfile_200bp > $bedfile_200bp'.output.bed'

### bigwig to tabs
###### covert reads count to NB p-value
while read LINE
do
	echo $LINE
	sig1=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $1}')
	sig2=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $2}')
	ct=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $3}')
	if ~/group/software/ucsc/bigWigAverageOverBed $bw_folder$sig1 $bedfile_200bp $bw_folder$sig1'.tab'; then echo 'bigWigAverageOverBed'; else echo 'ERROR: bigWigAverageOverBed' && exit 1; fi
	if sort -k1,1 $bw_folder$sig1'.tab' | cut -f5 > $bw_folder$sig1'.sig.tab'; then echo 'bigWigAverageOverBed'; else echo 'ERROR: bigWigAverageOverBed' && exit 1; fi
	if Rscript $script_dir'negative_binomial_p_2r_bgadj.R' $bw_folder$sig1'.sig.tab' 1 $bw_folder$ct'.'$sig1'.nbp.tab'; then echo 'covert reads count to NB p-value DONE'; else echo 'ERROR: covert reads count to NB p-value' && exit 1; fi
done < $bw_file_list

### fisher's method merge p-values
cut -f3 $bw_file_list | sort -u > ct_list.txt
for ct in $(cat ct_list.txt)
do
	echo $ct
	if Rscript $script_dir'fisher_pval.R' $ct '.nbp_2r_bgadj.txt' $bw_folder 100; then echo 'get Fisher method merged pval DONE'; else echo 'ERROR: get Fishers method merged pval' && exit 1; fi
done

### select reference sample
ls $bw_folder*.fisher_p.frip_snr.txt > fisherp.file_list.txt
if time Rscript $script_dir'get_mk_ref.R' fisherp.file_list.txt frip fisher.ref_frip.txt; then echo 'select reference dataset for s3norm'; else echo 'ERROR: select reference dataset for s3norm' && exit 1; fi

### normalize to all mk ref
### get 200bp win signal of allref
if ~/group/software/ucsc/bigWigAverageOverBed $bw_folder$bw_file_allref $bedfile_200bp $bw_folder$bw_file_allref'.tab'; then echo 'bigWigAverageOverBed allref'; else echo 'ERROR: bigWigAverageOverBed allref' && exit 1; fi
if sort -k1,1 $bw_folder$bw_file_allref'.tab' | cut -f5 > $bw_folder$bw_file_allref'.sig.tab'; then echo 'bigWigAverageOverBed allref'; else echo 'ERROR: bigWigAverageOverBed allref' && exit 1; fi

### normalize to all mk ref
mk_ref=$(head -1 fisher.ref_frip.txt)
if time python $script_dir's3norm_rotate_log_ref_mean.py' -n 500000 -a $bw_folder$bw_file_allref'.sig.tab' -b $mk_ref -u $upperlim1 -l $lowerlim1; then echo 's3norm normalize reference DONE'; else echo 'ERROR: s3norm normalize reference' && exit 1; fi


### normalize with mk
while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	sig2_celltype=$(echo "$LINE" | awk '{print $2}' | awk -F '.' -v OFS='\t' '{print $1"_"$2}')
	upperlim=100
	lowerlim=0
	echo $sig1 
	echo $sig2
	echo $sig2_celltype
	### set upper limit
	cat $sig1 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig1'.upperlim.txt'
	cat $sig2 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig2'.upperlim.txt' 
	### peak norm
	if time python $script_dir's3norm_rotate_log_z_mean.py' -n 500000 -a $sig1'.upperlim.txt' -b $sig2'.upperlim.txt' -u $upperlim -l $lowerlim; then echo 'pknorm across datasets DONE'; else echo 'ERROR: pknorm across datasets' && exit 1; fi
	### rm tmp files
	rm $sig1'.upperlim.txt'
	rm $sig2'.upperlim.txt'
done < fisher.ref_frip.txt.info.txt









###### convert nbp to Fisher's method merged p-value
### extrac cell type mark list
ls $bw_folder*'.nbp.tab' | awk -F '.' -v OFS='\t' '{print $1"."$2}' | sort -u > cell_marker_list.txt
ls *.nbp_2r_bgadj.txt | awk -F '.' -v OFS='\t' '{print $2}' | sort -u > mark_list.txt
ls *.nbp_2r_bgadj.txt | awk -F '.' -v OFS='\t' '{print $1}' | sort -u > cell_list.txt

### move data NB p-value data into nbp folder
if [ -d $working_dir'nbp/' ]; then echo $working_dir'nbp/' exist; else mkdir $working_dir'nbp/'; fi
mv *.nbp_2r_bgadj.txt $working_dir'nbp/'
mv *.mvsp.txt $working_dir'nbp/'
### get Fisher's method merged pval
for cm in $(cat cell_marker_list.txt)
do
	echo $cm
	if time Rscript $script_dir'fisher_pval.R' $cm '.nbp_2r_bgadj.txt' $working_dir'nbp/' 100; then echo 'get Fisher method merged pval DONE'; else echo 'ERROR: get Fishers method merged pval' && exit 1; fi
done



###### select reference dataset for pknorm
for mk in $(cat mark_list.txt)
do
	echo $mk
	ls *$mk*.frip_snr.txt > $mk'.file_list.txt'
	if time Rscript $script_dir'get_mk_ref.R' $mk'.file_list.txt' $select_method $mk'.ref_frip.txt'; then echo 'select reference dataset for pknorm'; else echo 'ERROR: select reference dataset for pknorm' && exit 1; fi
done
### select top reference dataset for cross mark pknorm
if time Rscript $script_dir'get_top_ref.R' '.ref_frip.txt' $select_method $working_dir cross_mark_ref_list.txt $user_given_global_ref; then echo 'select top reference dataset for cross mark pknorm DONE'; else echo 'ERROR: select top reference dataset for cross mark pknorm' && exit 1; fi



###### pknorm normalize reference datasets of all marks
while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	sig2_celltype=$(echo "$LINE" | awk '{print $2}' | awk -F '.' -v OFS='\t' '{print $1"_"$2}')
	upperlim=100
	lowerlim=0
	echo $sig1 
	echo $sig2
	echo $sig2_celltype
	### set upper limit
	cat $sig1 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig1'.upperlim.txt'
	cat $sig2 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig2'.upperlim.txt' 
	### peak norm
	if time python $script_dir's3norm_rotate_log_ref_mean.py' -n 500000 -a $sig1'.upperlim.txt' -b $sig2'.upperlim.txt' -u $upperlim -l $lowerlim; then echo 'pknorm normalize reference DONE'; else echo 'ERROR: pknorm normalize reference' && exit 1; fi
	### rm tmp files
	rm $sig1'.upperlim.txt'
	rm $sig2'.upperlim.txt'
done < cross_mark_ref_list.txt.info.txt
### move ref norm files into ref_info folder
if [ -d $working_dir'ref_info/' ]; then echo $working_dir'ref_info/' exist; else mkdir $working_dir'ref_info/'; fi
mv *.pknorm.scatterplot.png $working_dir'ref_info/'
mv *.scatterplot.png $working_dir'ref_info/'
mv *.ref.info.txt $working_dir'ref_info/'



###### pknorm across datasets with the same mark
for mk in $(cat mark_list.txt)
do
	echo $mk
	while read LINE
	do
		sig1=$(echo "$LINE" | awk '{print $1}')
		sig2=$(echo "$LINE" | awk '{print $2}')
		sig2_celltype=$(echo "$LINE" | awk '{print $2}' | awk -F '.' -v OFS='\t' '{print $1"_"$2}')
		upperlim=100
		lowerlim=0
		echo $sig1 
		echo $sig2
		echo $sig2_celltype
		### set upper limit
		cat $sig1 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig1'.upperlim.txt'
		cat $sig2 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig2'.upperlim.txt' 
		### peak norm
		if time python $script_dir's3norm_rotate_log_z_mean.py' -n 500000 -a $sig1'.upperlim.txt' -b $sig2'.upperlim.txt' -u $upperlim -l $lowerlim; then echo 'pknorm across datasets DONE'; else echo 'ERROR: pknorm across datasets' && exit 1; fi
		### rm tmp files
		rm $sig1'.upperlim.txt'
		rm $sig2'.upperlim.txt'
	done < $mk'.ref_frip.txt.info.txt'
done
### move pknorm files into pknorm_sig folder
if [ -d $working_dir'pknorm_info/' ]; then echo $working_dir'pknorm_info/' exist; else mkdir $working_dir'pknorm_info/'; fi
mv *.pknorm.scatterplot.png $working_dir'pknorm_info/'
mv *.scatterplot.png $working_dir'pknorm_info/'
mv *.info.txt $working_dir'pknorm_info/'



###### mv PKnorm normalized signal files & unnormalized signal files into *_sig folders
### pknorm signal
if [ -d $working_dir'pknorm_sig/' ]; then echo $working_dir'pknorm_sig/' exist; else mkdir $working_dir'pknorm_sig/'; fi
mv *.pknorm.txt $working_dir'pknorm_sig/'
### ref pknorm signal
if [ -d $working_dir'pknorm_ref_sig/' ]; then echo $working_dir'pknorm_ref_sig/' exist; else mkdir $working_dir'pknorm_ref_sig/'; fi
mv *.pknorm.ref.txt $working_dir'pknorm_ref_sig/'
### ref frip & snp
mv *.ref_frip.txt $working_dir'ref_info/'
### fisher pvalue signal without normalization
if [ -d $working_dir'fisherp/' ]; then echo $working_dir'fisherp/' exist; else mkdir $working_dir'fisherp/'; fi
mv *.fisher_p.txt $working_dir'fisherp/'
mv *.frip_snr.txt $working_dir'fisherp/'



###### set limit for signals
for filename in $(cat cell_marker_list.txt)
do
	echo $filename
	if cat $working_dir'pknorm_sig/'$filename'.pknorm.txt' | awk -F '\t' -v OFS='\t' -v ul=$overall_upper -v ll=$overall_lower '{if ($1<ll) print ll; else if ($1>ul) print ul; else print $1}' > $filename'.pknorm.'$overall_lower'_'$overall_upper'.txt'; then echo 'set limit for signals DONE'; else echo 'ERROR: set limit for signals' && exit 1; fi
done
### mv to the output folder
if [ -d $working_dir'pknorm_'$overall_lower'_'$overall_upper'_sig/' ]; then echo $working_dir'pknorm_'$overall_lower'_'$overall_upper'_sig/' exist; else mkdir $working_dir'pknorm_'$overall_lower'_'$overall_upper'_sig/'; fi
mv *'.pknorm.'$overall_lower'_'$overall_upper'.txt' $working_dir'pknorm_'$overall_lower'_'$overall_upper'_sig/'



### list files
if [ -d $working_dir'list_files/' ]; then echo $working_dir'list_files/' exist; else mkdir $working_dir'list_files/'; fi
mv *list.txt $working_dir'list_files/'





