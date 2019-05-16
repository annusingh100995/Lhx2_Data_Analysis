myarr2=($(ls -l|cut -c 46-))
  
 for j in "${myarr2[@]}"
	do
		echo ./$j 
		printf "In DIRECTORY $j \n"
		cd ./$j

		myarr=($(find ./ -name "*.fastq.gz"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
		printf "bla vla./${myarr[1]} ./${myarr[2]} \n"
		cd ../
	done


myarr5=($(ls -l|cut -c 46-))
  
 for j in "${myarr5[@]}"
	do
		printf "In DIRECTORY $j \n"
		cd ./$j
		fastqc *.fastq.gz 
		cd ../
	done


				if [ $a == *"HISAT2_INDEX"* ]
				then
	 			printf " HISAT_INDEX \n"
	 			continue   # break the for looop
				fi
 hisat2 -x /home/darwin/Lhx2_RNASeq/HISAT2_INDEX/hisat_indexed_file -1 ./E12.5_CTX_Lhx2_MUT_REP1/Lhx2KOE12-5ctx1_1.fastq.gz -2 ./E12.5_CTX_Lhx2_MUT_REP1/Lhx2KOE12-5ctx1_2.fastq.gz -S ./E12.5_CTX_Lhx2_MUT_REP1/E12.5_CTX_LHX2_REP1.sam -5 5 -3 10 -p 8 
#printf "bla vla./${myarr_replicates[1]} ./${myarr_replicates[2]} hisat_out_$j \n"

myarr_RNA=($(ls -l|cut -c 46-))
	for a in "${myarr_RNA[@]}"
			do
			cd ./$a
			printf "IN DIRECTORY JUST CHECKINMG $a \n " # either 12.5 or 15.5
			arr_stage=($(ls -l|cut -c 46-))
				for b in "${arr_stage[@]}"
					do 
					printf " IN DIRECTORY printing stage $a $b \n "  #either ctx or HC
					cd ./$b		
					myarr2_samples=($(ls -l|cut -c 46-))
  						for j in "${myarr2_samples[@]}"
							do
							echo ./$j 
							printf "In DIRECTORY SAMPLE NAMES  $j \n"
							cd ./$j
							myarr_replicates=($(find ./ -name "*.fastq.gz"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
	printf" hisat2  -p 8 -5 5 -3 10 -x /home/darwin/Lhx2_RNASeq/HISAT2_INDEX/hisat_indexed_file -1 ./${myarr_replicates[1]} -2 ./${myarr_replicates[2]} -S ./hisat_out_$j.sam" 
							cd ../
							done
					
					done
			
			done
							
	


#####################################################################################################################################################







MYARR_LHX2=($(ls -l|cut -c 46-))
for a in "${MYARR_LHX2[@]}"    
                cd ./$a # MOVES IN E12.5 /E15.5
                printf "IN DIRECTORY JUST CHECKINMG $a \n " # EITHER IN 12.5 or 15.5
                        MYARR_STAGE=($(ls -l|cut -c 46-))
                                for b in "${MYARR_STAGE${[@]}"
                                        do
                                        printf " IN DIRECTORY printing stage $a $b \n "  # EITHER CTX or HC
                                        cd ./$b # MOVES IN CTX / HC
                                        MYARR_SAMPLES=($(ls -l|cut -c 46-))
                                                for j in "${MYARR_SAMPLES[@]}"
                                                        do
                                                        echo ./$j
                                                        printf "In DIRECTORY SAMPLE NAMES  $j \n"
                                                        cd ./$j # MOVES IN ITH SAMPLE
                                                        MYARR_REPLICATES=($(find ./ -name "*.fastq.gz"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
                                                        printf "bla vla./${MYARR_REPLICATES[1]} ./${MYARR_REPLICATES[2]} \n"

                                                        cd ../ # MOVES OUT OF ITH SAMPLE
                                                        done
                                             cd ../ 
                                        done
				cd ../ 
                        done




######################################################################################################################################################################################


					 myarr1_samples_HC=($(ls -l|cut -c 46-))
                                                for j in "${myarr1_samples_HC[@]}"
                                                        do
                                                        printf "In DIRECTORY SAMPLE NAMES  $j \n"
                                                        cd ./$j
                                                        myarr_replicates_HC=($(find ./ -name "*.fastq.gz"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
                                                        printf "   ALIGNMENT OF Sample $j \n"
			 hisat2  -p 4 -5 5 -3 10 -x /home/darwin/Lhx2_RNASeq/HISAT2_INDEX/hisat_indexed_file -1 ./${myarr_replicates_HC[1]} -2 ./${myarr_replicates_HC[2]} -S ./hisat_out_$j.sam
                                                        cd ../-p
                                                        done







for a in "${MYARR_LHX2[@]}"   
	do 
        cd ./$a
	printf "IN DIRECTORY $a \n"
			MYARR_STAGE=($(ls -l|cut -c 46-))
			for b in "${MYARR_STAGE${[@]}"
                                        do
				cd ./$b
                                printf " IN DIRECTORY printing stage $a $b \n "  # EITHER CTX or HC
                                         
					
					cd ../
					done
		cd ../
done


MYARR_LHX2=($(ls -l|cut -c 46-))
for a in "${MYARR_LHX2[@]}"   
	do 
        cd ./$a
	printf "IN DIRECTORY $a \n"
			MYARR_STAGE=($(ls -l|cut -c 46-))
			for b in "${MYARR_STAGE[@]}"
			do
			cd ./$b
			printf " IN DIRECTORY printing stage $a  and $b \n "
			MYARR_SAMPLES=($(ls -l|cut -c 46-))
                                                for j in "${MYARR_SAMPLES[@]}"
                                                        do
                                                        echo ./$j
                                                        printf "In DIRECTORY SAMPLE NAMES  $j \n"
                                                        cd ./$j
                                                        MYARR_REPLICATES=($(find ./ -name "*.fastq.gz"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
                                                        printf "bla vla./${MYARR_REPLICATES[1]} ./${MYARR_REPLICATES[2]} \n"
							cd ../
							done

		
			cd ../
			done

	cd ../				
	done






featureCounts -p -T 8 -t gene -g gene_name -a /home/darwin/Lhx2_RNASeq/gene/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf -o ./LHX2_RNA_COUNT.txt E12.5/CTX/E12.5_CTX_Lhx2_MUT_REP1/hisat_out_E12.5_CTX_Lhx2_MUT_REP1.sam E12.5/CTX/E12.5_CTX_Lhx2_MUT_REP2/hisat_out_E12.5_CTX_Lhx2_MUT_REP2.sam E12.5/CTX/E12.5_CTX_Lhx2_MUT_REP3/hisat_out_E12.5_CTX_Lhx2_MUT_REP3.sam E12.5/CTX/E12.5_CTX_WT_REP1/hisat_out_E12.5_CTX_WT_REP1.sam E12.5/CTX/E12.5_CTX_WT_REP2/hisat_out_E12.5_CTX_WT_REP2.sam E12.5/CTX/E12.5_CTX_WT_REP3/hisat_out_E12.5_CTX_WT_REP3.sam E12.5/HC/E12.5_HC_LHX2_MUT_REP1/hisat_out_E12.5_HC_LHX2_MUT_REP1.sam E12.5/HC/E12.5_HC_LHX2_MUT_REP2/hisat_out_E12.5_HC_LHX2_MUT_REP2.sam E12.5/HC/E12.5_HC_LHX2_MUT_REP3/hisat_out_E12.5_HC_LHX2_MUT_REP3.sam E12.5/HC/E12.5_WT_HC_REP1/hisat_out_E12.5_WT_HC_REP1.sam E12.5/HC/E12.5_WT_HC_REP2/hisat_out_E12.5_WT_HC_REP2.sam E12.5/HC/E12.5_WT_HC_REP3/hisat_out_E12.5_WT_HC_REP3.sam E15.5/CTX/E15.5_CTX_LHX2_MUT_REP1/hisat_out_E15.5_CTX_LHX2_MUT_REP1.sam E15.5/CTX/E15.5_CTX_LHX2_MUT_REP2/hisat_out_E15.5_CTX_LHX2_MUT_REP2.sam E15.5/CTX/E15.5_CTX_LHX2_MUT_REP3/hisat_out_E15.5_CTX_LHX2_MUT_REP3.sam E15.5/CTX/E15.5_CTX_WT_REP1/hisat_out_E15.5_CTX_WT_REP1.sam E15.5/CTX/E15.5_CTX_WT_REP2/hisat_out_E15.5_CTX_WT_REP2.sam E15.5/CTX/E15.5_CTX_WT_REP3/hisat_out_E15.5_CTX_WT_REP3.sam E15.5/HC/E15.5_HC_LHX2_MUT_REP1/hisat_out_E15.5_HC_LHX2_MUT_REP1.sam E15.5/HC/E15.5_HC_LHX2_MUT_REP2/hisat_out_E15.5_HC_LHX2_MUT_REP2.sam E15.5/HC/E15.5_HC_LHX2_MUT_REP3/hisat_out_E15.5_HC_LHX2_MUT_REP3.sam E15.5/HC/E15.5_HC_WT_REP1/hisat_out_E15.5_HC_WT_REP1.sam E15.5/HC/E15.5_HC_WT_REP2/hisat_out_E15.5_HC_WT_REP2.sam E15.5/HC/E15.5_HC_WT_REP3/hisat_out_E15.5_HC_WT_REP3.sam










