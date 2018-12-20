# normalize_to_vision

### Normalize other signal to the s3norm normalized VISION signal

## clone the pipeline from github
```
git clone https://github.com/guanjue/normalize_to_vision.git
```

## requirements
###### python library
```
numpy
subprocess
matplotlib
scipy
```
###### R library
```
metap
LSD

```


## Input files
###### (1) the bw_file_list: (3 columns separated by tab)
###### 1st column: bw file name
###### 2nd column: always 1
###### 3rd column: cell type name. It will be used as the outputname (DO NOT put any '.' in this column!)
```
>head bw_file_list.txt
ideasVisionV20p8NormCd4Atac.bw	1	CD4
ideasVisionV20p8NormCd8Atac.bw	1	CD8
ideasVisionV20p8NormCmpAtac.bw	1	CMP
ideasVisionV20p8NormCmpAtac1.bw	1	CMP
```

###### (2) The bed file of 200bp window with 4 columns
###### !!! remove blacklist first !!!
```
head mm10.200bin.rand.sort.bed
chr1	2274200	2274400	chr1_2274200_2274400
chr1	5290600	5290800	chr1_5290600_5290800
chr1	6490400	6490600	chr1_6490400_6490600
chr1	6645600	6645800	chr1_6645600_6645800
chr1	6718000	6718200	chr1_6718000_6718200
chr1	7196000	7196200	chr1_7196000_7196200
chr1	7837000	7837200	chr1_7837000_7837200
chr1	8444000	8444200	chr1_8444000_8444200
chr1	8789200	8789400	chr1_8789200_8789400
chr1	8795200	8795400	chr1_8795200_8795400
```

###### (3) All of the bw file needs to be save in the $bw_folder folder

### !!! (4) !!!
###### (4) overall ref bw file: The signal track used to scale different marks
###### Always use the ideasVisionV20p8NormImkK27ac.bw as the reference between marks
```
wget http://usevision.org/data/IDEASmouseHem/PKnorm/ideasVisionV20p8NormImkK27ac.bw 
mv ideasVisionV20p8NormImkK27ac.bw user_given_bw_folder/
```


## RUN pipeline
##### (1) Copy run_normalize_to_vision.sh to the working directory. Change the folder names (script_dir, bw_folder, output_folder, bedfile_200bp_filename_remove_blacklist, bw_file_list) in 'run_normalize_to_vision.sh'
##### The folder path should always be the absolute path
```
##################################
script_dir='/storage/home/g/gzx103/group/projects/vision/normalize_to_vision/src/'
bw_folder='/storage/home/g/gzx103/group/projects/vision/normalize_to_vision/signals/'
output_folder='/storage/home/g/gzx103/group/projects/vision/normalize_to_vision/s3norm_output/'
bedfile_200bp='mm10.200bin.rand.sort.bed'
bw_file_list='bw_file_list.txt'
```

##### (2) create the 'bw_file_list.txt' file in the working directory
```
>head bw_file_list.txt
ideasVisionV20p8NormCd4Atac.bw	1	CD4
ideasVisionV20p8NormCd8Atac.bw	1	CD8
ideasVisionV20p8NormCmpAtac.bw	1	CMP
ideasVisionV20p8NormCmpAtac1.bw	1	CMP
```

##### (3) use 'run_normalize_to_vision.sh' script to run pipeline
```
time bash run_normalize_to_vision.sh
```

## Output files
##### The output files will be saved in the 'output_folder'
##### (1) The bed file of all 200bp-bins: 's3norm_output/'$bedfile_200bp'.output.bed'
##### The peak is sorted by peak name (4th column of the file)
##### For running IDEAS, the signal of other marks of this 200-bins bed file should also be generated

##### (2) The s3norm normalized signal of each cell type of mark named as $ct'.s3norm.txt'
##### The each file only has one column that has the s3norm normalized signal of each 200bp-bin.
##### The row order is the same as the 's3norm_output/'$bedfile_200bp'.output.bed' file
```
>ls -ltrh s3norm_output/*.s3norm.txt
-rw-r-----  1 universe  staff    97K Dec 20 17:40 s3norm_output/CMP.s3norm.txt
-rw-r-----  1 universe  staff    98K Dec 20 17:40 s3norm_output/CD4.s3norm.txt
>
>
>
>head s3norm_output/CMP.s3norm.txt
0.04674804690559381
0.04674804690559381
0.04674804690559381
0.04674804690559381
```







