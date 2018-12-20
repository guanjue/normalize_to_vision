###### read parameters from inputs
script_dir=$1
bw_folder=$2
output_folder=$3
bedfile_200bp=$4
bw_file_list=$5
bw_file_allref=$6
upperlim1=$7
lowerlim1=$8

### get testing data
#wget http://usevision.org/data/IDEASmouseHem/PKnorm/ideasVisionV20p8NormCmpAtac.bw
#wget http://usevision.org/data/IDEASmouseHem/PKnorm/ideasVisionV20p8NormCd4Atac.bw
#wget http://usevision.org/data/IDEASmouseHem/PKnorm/ideasVisionV20p8NormCd8Atac.bw
#cp ideasVisionV20p8NormCmpAtac.bw ideasVisionV20p8NormCmpAtac1.bw ### as a replicate of CMPatac

### get random windows
#bedtools makewindows -g ~/group/genome/mm10/mm10.chrom.sizes -w 200 > mm10.200bin.bed
#Rscript get_random_200bp_bins.R
#cat mm10.200bin.rand.bed | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' | sort -k1,1 -k2,2n > $bedfile_200bp

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
	if time python $script_dir's3norm_rotate_log_z_mean.py' -n 500000 -a $sig1'.upperlim.txt' -b $sig2'.upperlim.txt' -u $upperlim -l $lowerlim; then echo 's3norm across datasets DONE'; else echo 'ERROR: s3norm across datasets' && exit 1; fi
	### rm tmp files
	rm $sig1'.upperlim.txt'
	rm $sig2'.upperlim.txt'
done < fisher.ref_frip.txt.info.txt


### mv all signal files to output folder
if [ -d $output_folder ]; then echo $output_folder exist; else mkdir $output_folder; fi
mv $bw_folder*'.s3norm.txt' $output_folder
mv $bedfile_200bp'.output.bed' $output_folder
mv ct_list.txt $output_folder
mv fisherp.file_list.txt $output_folder
mv fisher.ref_frip.txt $output_folder
mv fisher.ref_frip.txt.info.txt $output_folder

