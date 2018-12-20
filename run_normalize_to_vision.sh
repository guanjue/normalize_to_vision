###### set parameters
script_dir='/storage/home/g/gzx103/group/projects/vision/normalize_to_vision/src/'
bw_folder='/storage/home/g/gzx103/group/projects/vision/normalize_to_vision/signals/'
output_folder='/storage/home/g/gzx103/group/projects/vision/normalize_to_vision/s3norm_output/'
bedfile_200bp='mm10.200bin.rand.sort.bed'
bw_file_list='bw_file_list.txt'
bw_file_allref='ideasVisionV20p8NormImkK27ac.bw'
upperlim1=16
lowerlim1=0

time bash $script_dir'normailize_to_vision_overall_pipeline.sh' $script_dir $bw_folder $output_folder $bedfile_200bp $bw_file_list $bw_file_allref $upperlim1 $lowerlim1

