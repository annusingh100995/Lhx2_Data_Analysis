

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
cuffmerge -g /home/darwin/Lhx2_RNASeq/gene/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf -s /home/darwin/Lhx2_RNASeq/genome/GRCm38.p5.genome.fa -p 4 -o ./merged_files/mergerd_$p ./$j.txt
done

##########################CUFFDIFF########################################################


cuffdiff -o diff_out_sampledata -b ../genome/GRCm38.p5.genome.fa -p 1 -L CONTROL,RAS -u ./merged_asm/merged.gtf ./CONTROL_REP1_thout/accepted_hits.bam,CONTROL_REP2_thout/accepted_hits.bam,CONTROL_REP3_thout/accepted_hits.bam ./RAS_REP1_thout/accepted_hits.bam,RAS_REP2_thout/accepted_hits.bam,RAS_REP3_thout/accepted_hits.bam 

cuffdiff -o ./diff_out/diff_out_E12.5_E15.5_WT_ctx -b ./genome/GRCm38.p5.genome.fa -p 4 -L E12.5_WT_ctx,E15.5_WT_ctx -u ./assemblies/merged_files/mergerd_amb_E12.5_E15.5_WT_ctx/merged.gtf ./E12.5/E12.5_cortex/thout_E12.5_ctx/thout_E12.5_WT_ctx_rep1_CGATGT/accepted_hits.bam,./E12.5/E12.5_cortex/thout_E12.5_ctx/thout_E12.5_WT_ctx_rep2_TTAGGC/accepted_hits.bam ./E15.5/E15.5_cortex/thout_E15.5_ctx/thout_E15.5_WT_ctx_rep1/accepted_hits.bam,./E15.5/E15.5_cortex/thout_E15.5_ctx/thout_E15.5_WT_ctx_rep2/accepted_hits.bam 

cuffdiff -o ./diff_out/diff_out_E12.5_E15.5_WT_HC -b ./genome/GRCm38.p5.genome.fa -p 4 -L E12.5_WT_HC,E15.5_WT_HC -u ./assemblies/merged_files/mergerd_amb_E12.5_E15.5_WT_HC/merged.gtf ./E12.5/E12.5_HC/thout_E12.5_HC/thout_E12.5_WT_HC_rep1_ACTTGA/accepted_hits.bam,./E12.5/E12.5_HC/thout_E12.5_HC/thout_E12.5_WT_HC_rep1_CATGATC/accepted_hits.bam ./E15.5/E15.5_HC/thout_E15.5_HC/thout_E15.5_WT_HC_rep1/accepted_hits.bam,./E15.5/E15.5_HC/thout_E15.5_HC/thout_E15.5_WT_HC_rep2/accepted_hits.bam 

cuffdiff -o ./diff_out/diff_out_E12.5_WT_ctx_HC -b ./genome/GRCm38.p5.genome.fa -p 4 -L E12.5_WT_ctx,12.5_WT_HC ./assemblies/merged_files/mergerd_amb_E12.5_WT_ctx_HC/merged.gtf ./E12.5/E12.5_cortex/thout_E12.5_ctx/thout_E12.5_WT_ctx_rep1_CGATGT/accepted_hits.bam,./E12.5/E12.5_cortex/thout_E12.5_ctx/thout_E12.5_WT_ctx_rep2_TTAGGC/accepted_hits.bam ./E12.5/E12.5_HC/thout_E12.5_HC/thout_E12.5_WT_HC_rep1_ACTTGA/accepted_hits.bam,./E12.5/E12.5_HC/thout_E12.5_HC/thout_E12.5_WT_HC_rep1_CATGATC/accepted_hits.bam

cuffdiff -o ./diff_out/diff_out_E12.5_WT_Lhx2_mut_ctx -b ./genome/GRCm38.p5.genome.fa -p 4 -L E12.5_WT_ctx,E12.5_Lhx2_mut_ctx ./assemblies/merged_files/merged_amb_E12.5_WT_Lhx2_mut_ctx/merged.gtf ./E12.5/E12.5_cortex/thout_E12.5_ctx/thout_E12.5_WT_ctx_rep1_CGATGT/accepted_hits.bam,./E12.5/E12.5_cortex/thout_E12.5_ctx/thout_E12.5_WT_ctx_rep2_TTAGGC/accepted_hits.bam ./E12.5/E12.5_cortex/thout_E12.5_ctx/thout_E12.5_Lhx2_mut_ctx_rep1_ACAGTG/accepted_hits.bam,./E12.5/E12.5_cortex/thout_E12.5_ctx/thout_E12.5_Lhx2_mut_ctx_rep2_GCCAAT/accepted_hits.bam

cuffdiff -o ./diff_out/diff_out_E12.5_WT_lhx2_mut_HC -b ./genome/GRCm38.p5.genome.fa -p 4 -L E12.5_WT_HC,E12.5_Lhx2_mut_HC ./assemblies/merged_files/mergerd_amb_E12.5_WT_Lhx2_mut_HC/merged.gtf ./E12.5/E12.5_HC/thout_E12.5_HC/thout_E12.5_WT_HC_rep1_ACTTGA/accepted_hits.bam,./E12.5/E12.5_HC/thout_E12.5_HC/thout_E12.5_WT_HC_rep1_CATGATC/accepted_hits.bam ./E12.5/E12.5_HC/thout_E12.5_HC/thout_E12.5_Lhx2_mut_HC_rep1_GGCTAC/accepted_hits.bam,./E12.5/E12.5_HC/thout_E12.5_HC/thout_E12.5_Lhx2_mut_HC_rep2_CTTGTA/accepted_hits.bam

cuffdiff -o ./diff_out/diff_out_E15.5_WT_ctx_HC -b ./genome/GRCm38.p5.genome.fa -p 4 -L E15.5_WT_ctx,E15.5_WT_HC ./assemblies/merged_files/mergerd_amb_E15.5_WT_ctx_HC/merged.gtf ./E15.5/E15.5_cortex/thout_E15.5_ctx/thout_E15.5_WT_ctx_rep1/accepted_hits.bam,./E15.5/E15.5_cortex/thout_E15.5_ctx/thout_E15.5_WT_ctx_rep2/accepted_hits.bam ./E15.5/E15.5_HC/thout_E15.5_HC/thout_E15.5_WT_HC_rep1/accepted_hits.bam,./E15.5/E15.5_HC/thout_E15.5_HC/thout_E15.5_WT_HC_rep2/accepted_hits.bam 

cuffdiff -o diff_out/diff_out_E15.5_WT_Lhx2_mut_ctx -b ./genome/GRCm38.p5.genome.fa -p 4 -L E15.5_WT_ctx,E15.5_Lhx2_mut_ctx ./assemblies/merged_files/mergerd_amb_E15.5_WT_Lhx2_mut_ctx/merged.gtf ./E15.5/E15.5_cortex/thout_E15.5_ctx/thout_E15.5_WT_ctx_rep1/accepted_hits.bam,./E15.5/E15.5_cortex/thout_E15.5_ctx/thout_E15.5_WT_ctx_rep2/accepted_hits.bam ./E15.5/E15.5_cortex/thout_E15.5_ctx/thout_E15.5_Lhx2_mut_ctx_rep1/accepted_hits.bam,./E15.5/E15.5_cortex/thout_E15.5_ctx/thout_E15.5_Lhx2_mut_ctx_rep2/accepted_hits.bam 
 
cuffdiff -o diff_out/diff_out_E15.5_WT_Lhx2_mut_HC -b ./genome/GRCm38.p5.genome.fa -p 4 -L E15.5_WT_HC,E15.5_Lhx2_mut_HC ./assemblies/merged_files/mergerd_amb_E15.5_WT_Lhx2_mut_HC/merged.gtf ./E15.5/E15.5_HC/thout_E15.5_HC/thout_E15.5_WT_HC_rep1/accepted_hits.bam,./E15.5/E15.5_HC/thout_E15.5_HC/thout_E15.5_WT_HC_rep2/accepted_hits.bam ./E15.5/E15.5_HC/thout_E15.5_HC/thout_E15.5_Lhx2_mut_HC_rep1/accepted_hits.bam,./E15.5/E15.5_HC/thout_E15.5_HC/thout_E15.5_Lhx2_mut_HC_rep2/accepted_hits.bam 

featureCounts -p -T 1 -a gene/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf -t gene -g gene_name -o COUNTS/E15.5_WT_Lhx2_mut_HC.txt ./E15.5/E15.5_HC/thout_E15.5_HC/thout_E15.5_WT_HC_rep1/accepted_hits.bam ./E15.5/E15.5_HC/thout_E15.5_HC/thout_E15.5_WT_HC_rep2/accepted_hits.bam ./E15.5/E15.5_HC/thout_E15.5_HC/thout_E15.5_Lhx2_mut_HC_rep1/accepted_hits.bam ./E15.5/E15.5_HC/thout_E15.5_HC/thout_E15.5_Lhx2_mut_HC_rep2/accepted_hits.bam 

featureCounts -p -T 1 -a gene/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf -t gene -g gene_name -o COUNTS/E15.5_WT_Lhx2_mut_ctx.txt ./E15.5/E15.5_cortex/thout_E15.5_ctx/thout_E15.5_WT_ctx_rep1/accepted_hits.bam ./E15.5/E15.5_cortex/thout_E15.5_ctx/thout_E15.5_WT_ctx_rep2/accepted_hits.bam ./E15.5/E15.5_cortex/thout_E15.5_ctx/thout_E15.5_Lhx2_mut_ctx_rep1/accepted_hits.bam ./E15.5/E15.5_cortex/thout_E15.5_ctx/thout_E15.5_Lhx2_mut_ctx_rep2/accepted_hits.bam 

 
