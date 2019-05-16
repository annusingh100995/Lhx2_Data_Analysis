

#!/bin/bash

#############TOPHAT############
myarr2=($(ls -l|cut -c 46-))
  
 for j in "${myarr2[@]}"
	do
		echo ./$j 
		cd ./$j
		ls
		
		myarr=($(find ./ -name "*.fq"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
			echo ${myarr[1]} ${myarr[2]}
		printf "tophat2 -p 4 /home/darwin/Lhx2_RNASeq/gene/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf -o /home/darwin/Lhx2_RNASeq/E12.5/E12.5_cortex/E12.5_Lhx2_mut_ctx_rep1_ACAGTG_thout /home/darwin/Lhx2_RNASeq/genome/mm10_indexed_file ./${myarr[1]} ./${myarr[2]}" 


tophat -p 4 -G /home/darwin/Lhx2_RNASeq/gene/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf -o /home/darwin/Lhx2_RNASeq/E12.5/E12.5_cortex/thout/$i /home/darwin/Lhx2_RNASeq/genome/mm10_indexed_file /home/darwin/


		cd ../

###############CUFFLINKS##################3




myarr2=($(ls -l|cut -c 46-))
  
 for j in "${myarr2[@]}"
	do
		echo ./$j 
		cd ./$j
		ls
		myarr=($(find ./ -name "*.bam"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
			echo ${myarr[1]}  ${myarr[2]}
		printf "IN DIRECTORY $j"
		cufflinks -p 4 -o clout_$j ./${myarr[2]}
		cd ../
	done

######################RENAME STUPID FILE NAME ERROR###########################3
myarr2=($(ls -l|cut -c 46-))
  
 for j in "${myarr2[@]}"
	do
		echo ./$j 
		cd ./$j
		prename 's/clout_thout/clout/' clout_thout*/
		cd ../
	done

########################ASSEMBLIES AND CUFFMERGE######################################
myarr2=($(find ./ -name "*.txt"|awk 'BEGIN {FS ="/"}{print$NF}'|rev| cut -c 5-|rev))|
for j in "${myarr2[@]}"
	do
p=$j|rev| cut -c 5-|rev
printf "cuffmerge -g ./gene/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf -s ./genome/GRCm38.p5.genome.fa -p 4 -o ./merged_files/mergerd_$p ./assemblies/$j.txt"
done









 
