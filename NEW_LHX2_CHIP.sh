
				myarr1=($(ls -l|cut -c 46-))
                                                for j in "${myarr1[@]}"
                                                        do
                                                        printf "In DIRECTORY SAMPLE NAMES  $j \n"
                                                        cd ./$j
                                                        myarr_replicates=($(find ./ -name "*.fastq.gz"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
                                                        printf "   ALIGNMENT OF Sample $j \n"
	 bowtie2 -p 4 --very-sensitive-local -X 300 -x /home/darwin/Lhx2_RNASeq/genome/mm10_indexed_file -1 ./${myarr_replicates[1]} -2 ./${myarr_replicates[2]} -S ./bowtie_out_$j.sam
			 				cd ../
                                                        done











hisat2  -p 4 -5 5 -3 10 -x /home/darwin/Lhx2_RNASeq/HISAT2_INDEX/hisat_indexed_file -1 ./${myarr_replicates[1]} -2 ./${myarr_replicates[2]} -S ./hisat_out_$j.sam
                                                        printf "bla vla./${MYARR_REPLICATES[1]} ./${MYARR_REPLICATES[2]} \n"


MYARR_STAGE_e15_5=($(ls -l|cut -c 46-))
			for b in "${MYARR_STAGE_e15_5[@]}"
			do
			cd ./$b
			printf " IN DIRECTORY printing stage $b \n "
			MYARR_SAMPLES_e15_5=($(ls -l|cut -c 46-))
                                                for j in "${MYARR_SAMPLES_e15_5[@]}"
                                                        do
                                                        echo ./$j
                                                        printf "In DIRECTORY SAMPLE NAMES  $j \n"
                                                        cd ./$j
                                                        MYARR_REPLICATES_e15_5=($(find ./ -name "*.fastq.gz"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
bowtie2 -p 4 --very-sensitive-local -X 300 -x /home/darwin/Lhx2_RNASeq/genome/mm10_indexed_file -1 ./${MYARR_REPLICATES_e15_5[1]} -2 ./${MYARR_REPLICATES_e15_5[2]} -S ./bowtie_out_$j.sam

							cd ../
							done

		
			cd ../
			done




MYARR_foldername=($(ls -l|cut -c 51-|rev|cut -c 12-|rev))
for i in "${MYARR_foldername[@]}"
	do 
	mkdir ./$i
	done


MYARR_HISTONE=($(ls -l|cut -c 46-))
 for i in "${MYARR_HISTONE[@]}"
	do 
	cd ./$i 
	ls -l 
 		MYARR_FOLDER2=($(ls -l|cut -c 46-))
			for j in "${MYARR_FOLDER2[@]}"
			do 	
			cd ./$j
			ls -l 
				MYARR_Folder3=($(ls -l|cut -c 46-))
					for l in "${MYARR_Folder3[@]}"
						do 
						cd ./$l
						MYARR_SAMPLES=($(find ./ -name "*.fastq.gz"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
						ls -l 			
						printf "cat_$j"
						cd ../
						done
			cd ../
			done
	cd ../
	done


MYARR_HISTONE=($(ls -l|cut -c 46-))
 for i in "${MYARR_HISTONE[@]}"
	do 
	cd ./$i 
MYARR_SAMPLES=($(find ./ -name "*.fastq.gz"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
			printf " cat ./${MYARR_SAMPLES[1]} ./${MYARR_SAMPLES[2]} > ./cat_$i.fastq \n"
						
			 cat ./${MYARR_SAMPLES[1]} ./${MYARR_SAMPLES[2]} > ./cat_$i.fastq 
						cd ../
						done
			

cat ./${MYARR_SAMPLES_IP[1]} ./${MYARR_SAMPLES_IP[2]} > ./cat_$i 
MYARR_Folder3=($(ls -l|cut -c 46-))
				for l in "${MYARR_Folder3[@]}"
				do 
				cd ./$l

$ printf '%s\n' "${my_array[@]}" 
cat ./${MYARR_SAMPLES_IP[1]} ./${MYARR_SAMPLES_IP[2]} > ./cat_$j.fastq

MYARR_HISTONE_IP=($(ls -l|cut -c 46-))
 for i in "${MYARR_HISTONE_IP[@]}"
	do 
	cd ./$i 
printf " DIRECTORY: $i \n"
MYARR_FOLDER2_IP=($(ls -l|cut -c 46-))
			for j in "${MYARR_FOLDER2_IP[@]}"
			do 	
			cd ./$j
printf "DIRECTORY: $j \n"
MYARR_SAMPLES_IP=($(find ./ -name "*.fastq.gz"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
			printf " cat ./${MYARR_SAMPLES_IP[1]} ./${MYARR_SAMPLES_IP[2]} > ./cat_$j.fastq \n"
			cat ./${MYARR_SAMPLES_IP[1]} ./${MYARR_SAMPLES_IP[2]} > ./cat_$j.fastq			
			 
			cd ../
			done

cd ../
done




samples=($(ls -l|cut -c 46-))
for i in "${samples[@]}"
	do 
	cd ./$i 
	con_samples=($(find ./ -name "*.fastq.gz"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
		printf " DIRECTORY : $i \n "
			printf " ./${con_samples[1]} ./${con_samples[2]} \n"
 bowtie2 -p 4 --very-sensitive-local -X 300 -x /home/darwin/Lhx2_RNASeq/genome/mm10_indexed_file -1 ./${con_samples[1]} -2 ./${con_samples[2]} -S ./bowtie_out_$i.sam

	cd ../
	done





samples=($(ls -l|cut -c 46-))
for i in "${samples[@]}"
	do 
	cd ./$i 
	con_samples=($(find ./ -name "*.sam"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
		printf " ./${con_samples[1]} \n "
	 printf  " samtools view -S -b ./${con_samples[1]} > ./$i.bam \n "
		samtools view -S -b ./${con_samples[1]} > ./$i.bam 
	cd ../
	done











MYARR_STAGE_e15_5=($(ls -l|cut -c 46-))
			for b in "${MYARR_STAGE_e15_5[@]}"
			do
			cd ./$b
			printf " IN DIRECTORY printing stage $b \n "
			MYARR_SAMPLES_e15_5=($(ls -l|cut -c 46-))
                                                for j in "${MYARR_SAMPLES_e15_5[@]}"
                                                        do
                                                        echo ./$j
                                                        printf " In DIRECTORY SAMPLE NAMES  $j \n"
                                                        cd ./$j
                                                        MYARR_REPLICATES_e15_5=($(find ./ -name "*.sam"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
							printf  " samtools view -S -b ./${MYARR_REPLICATES_e15_5[1]} > ./$j.bam \n "
							samtools view -S -b ./${MYARR_REPLICATES_e15_5[1]} > ./$i.bam 
							cd ../
							done

		
			cd ../
			done









samtools view -S -b ./${MYARR_SAMPLES_IP[1]} > ./$i.bam 
				samtools view -S -b ./${MYARR_SAMPLES_IP[1]} > ./$a.bam

MYARR_HISTONE_IP=($(ls -l|cut -c 46-))
 for a in "${MYARR_HISTONE_IP[@]}"
	do 
	cd ./$a
printf " DIRECTORY: $a \n"

MYARR_SAMPLES_IP=($(find ./ -name "*.sam"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
			
			printf  " samtools view -S -b ./${MYARR_SAMPLES_IP[1]} > ./$a.bam \n" 
			samtools view -S -b ./${MYARR_SAMPLES_IP[1]} > ./$a.bam
			cd ../
			done






find / -type f -exec grep -H 'text-to-find-here' {} \;


MYARR_STAGE_e12_5=($(ls -l|cut -c 46-))




macs14 -t E10.5_LHX2_IP_REP1/E10.5_LHX2_IP_REP1.bam -c E10.5_INPUT_REP1/E10.5_INPUT_REP1.bam -f BAM -g mm -n ./E10.5_INPUT_LHX2_IP_REP1
macs14 -t E10.5_LHX2_IP_REP2/E10.5_LHX2_IP_REP2.bam -c E10.5_INPUT_REP2/E10.5_INPUT_REP2.bam -f BAM -g mm -n ./E10.5_INPUT_LHX2_IP_REP2











