macs14 -t E10.5_LHX2_IP_REP1/E10.5_LHX2_IP_REP1.bam -c E10.5_INPUT_REP1/E10.5_INPUT_REP1.bam -f BAM -g mm -n ./E10.5_INPUT_LHX2_IP_REP1
macs14 -t E10.5_LHX2_IP_REP2/E10.5_LHX2_IP_REP2.bam -c E10.5_INPUT_REP2/E10.5_INPUT_REP2.bam -f BAM -g mm -n ./E10.5_INPUT_LHX2_IP_REP2



 samtools view -S -b ./bowtie_out_E15.5_CTX_LHX2_MUT_REP1.sam > ./E15.5_CTX_LHX2_REP1.bam 
 samtools view -S -b ./bowtie_out_E15.5_CTX_LHX2_MUT_REP2.sam > ./E15.5_CTX_LHX2_REP2.bam 
 samtools view -S -b ./bowtie_out_E15.5_HC_INPUT_REP1.sam > ./E15.5_HC_INPUT_REP1.bam 
 samtools view -S -b ./bowtie_out_E15.5_HC_INPUT_REP2.sam > ./E15.5_HC_INPUT_REP2.bam 

 samtools view -S -b ./bowtie_out_E15.5_HC_LHX2_MUT_REP1.sam > ./E15.5_HC_LHX2_REP1.bam 

 samtools view -S -b ./bowtie_out_E15.5_HC_LHX2_MUT_REP2.sam > ./E15.5_HC_LHX2_REP2.bam 



macs14 -t E10.5_LHX2_IP_REP2/E10.5_LHX2_IP_REP2.bam -c E10.5_INPUT_REP2/E10.5_INPUT_REP2.bam -f BAM -g mm -n E10.5_INPUT_LHX2_MUT_REP2 --single-profile --wig --space=25

macs14 -t E12.5_CTX_LHX2_REP1/E12.5_CTX_LHX2_REP1.bam -c E12.5_CTX_INPUT_REP1/E12.5_CTX_INPUT_REP1.bam -f BAM -g mm -n E12.5_CTX_INPUT_LHX2_REP1 --single-profile --wig --space=25




myarr2=($(ls -l|cut -c 46-))
  
 for j in "${myarr2[@]}"
	do
		echo ./$j 
		printf "In DIRECTORY $j \n"
		cd ./$j

		myarr=($(find ./ -name "*.fastq.gz"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
		printf "ALIGNING SAMPLES ./${myarr[1]} ./${myarr[2]} \n"
	printf "bowtie2 -p 8 --very-sensitive-local -X 300 -x /home/darwin/Lhx2_RNASeq/genome/mm10_indexed_file -1 ./${myarr[1]} -2 ./${myarr[2]} -S $j.sam \n "
	bowtie2 -p 8 --very-sensitive-local -X 300 -x /home/darwin/Lhx2_RNASeq/genome/mm10_indexed_file -1 ./${myarr[1]} -2 ./${myarr[2]} -S $j.sam		
	cd ../
	done


myarr2=($(ls -l|cut -c 46-))
  
 for j in "${myarr2[@]}"
	do 
		printf "In DIRECTORY $j \n"
		cd ./$j

		myarr=($(find ./ -name "*.sam"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
		rm ./${myarr[1]}
		cd ..
		done






SORTING AND INDENING

samtools sort ./${myarr[1]} ./ $j.sorted 
samtools index $j.sorted.bam 


myarr2=($(ls -l|cut -c 46-))
  
 for j in "${myarr2[@]}"
	do
		echo ./$j 
		printf "In DIRECTORY $j \n"
		cd ./$j

		myarrs=($(find ./ -name "*.bam"|awk 'BEGIN {FS ="/"}{print$NF}'| cut -c 1-))
		printf " samtools sort ./${myarrs[1]} $j.sorted \n grep .bam \n "		
		                                                       
		printf " samtools index $j.sorted.bam \n"
		
		cd ../
		done


	samtools sort $i $flag.sorted
myarr2=($(find ./ -name "*.bam"))
for i in "${myarr2[@]}"
 do
        flag=$( echo "$i" | rev |cut -c 5-|rev ) 
        printf " samtools sort $i $flag.sorted.bam \n "
	samtools sort $i $flag.sorted.bam
done


myarr3=($(find ./ -name "*.sorted.bam"))
for j in "${myarr3[@]}"
	do 
	samtools index $j 
	done



MYARR12=($(find ./ -name "*.sam"))
for i in "${MYARR12[@]}"
	do 
	flag=$( echo "$i" | rev |cut -c 5-21|rev )
	printf " flag \n"
	printf "$i \n"
	printf "samtools -S -b $i > $flag.bam \n"
	done

 

sed -n -e 's/^.*'

samfiles=($(find ./ -name "*.sam"))    
for i in "${samfiles[@]}"
do 
	flag=$(echo "$i"|rev|cut -c 5-|rev)
	printf " $flag \n "
done




samfiles=($(find ./ -name "*.sam"))    
for i in "${samfiles[@]}"
do 
	flag=$(echo "$i"|rev|cut -c 5-|rev| cut -c 7-|sed 's/*//')
	printf	" $flag \n"
done



myarr2=($(find ./ -name "*.sam"))



for i in "${myarr2[@]}"
 do
        flag=$( echo "$i" | rev |cut -c 5-|rev ) 
        printf " samtools view -S -b $i > $flag.bam \n "
	samtools view -S -b $i > $flag.bam 
done



myarr_1=($(find ./ -name "*.fastq.gz"))
for i in "${myarr_1[@]}"
do
	flag=$( echo "$i" | rev |cut -c 12-|rev|cut -c 3-)
	printf " $flag \n"
	mkdir $flag
done 


myarr27=($(find ./ -name "*.sorted.bam"))
for i in "${myarr2[@]}"
 do
	flag=$( echo "$i" | rev |cut -c 12-|rev)

	printf "bamCoverage --normalizeUsingRPKM --extendReads 200 -b $i -o $i.bw \n "
	 bamCoverage --normalizeUsingRPKM --extendReads 200 -b $i -o $i.bw
done






myarr2=($(find ./ -name "*.sorted.bam"))
for j in "${myarr2[@]}"
	do 
	printf "samtools index $j \n "
	samtools index $j
	done


myarr_input=($(find ./ -name "*INPUT*.sorted.bam"))

myarr_Lhx2=($(find ./ -name "*LHX2*.sorted.bam"))

for ((i=1;i<=5;i++)); do
        printf " ${myarr_input[$i]} : ${myarr_Lhx2[$i]} \n "
done



HISTONE
 

myarr_input=($(find ./ -name "*cnt*.sorted.bam"))

myarr_Lhx2=($(find ./ -name "*lhx2-mut*.sorted.bam"))

for ((i=1;i<=4;i++)); do
	flag=$( echo "${myarr_Lhx2[$i]}" | rev |cut -c 12-|rev)
	printf " plotFingerprint -b ${myarr_Lhx2[$i]} ${myarr_input[$i]} -plot $flag.png \n"        
	plotFingerprint -b ${myarr_Lhx2[$i]} ${myarr_input[$i]} -plot $flag.png
done

CTX

LHX2_HISTONE_CTX=($(find ./ -name "*LHX2_MUT*H*.sorted.bam"))
for j in "${LHX2_HISTONE_CTX[@]}"
	do 
	flag=$( echo "$j" | rev |cut -c 12-|rev)
	printf " plotFingerprint -b $j ./E12.5-CTX-INPUT/E12.5_CORTEX_INPUT_LHX2_MUT/E12.5_CORTEX_INPUT_LHX2_MUT.sorted.bam -plot $flag.png \n"
	plotFingerprint -b $j ./E12.5-CTX-INPUT/E12.5_CORTEX_INPUT_LHX2_MUT/E12.5_CORTEX_INPUT_LHX2_MUT.sorted.bam -plot $flag.png	
	done



CNT_HISTONE_CTX=($(find ./ -name "*CNT*H*.sorted.bam"))
for j in "${CNT_HISTONE_CTX[@]}"
	do 
	flag=$( echo "$j" | rev |cut -c 12-|rev)
	printf " plotFingerprint -b $j  E12.5-CTX-INPUT/E12.5_CORTEX_INPUT_CNT/E12.5_CORTEX_INPUT.sorted.bam -plot $flag.png \n"
	plotFingerprint -b $j  E12.5-CTX-INPUT/E12.5_CORTEX_INPUT_CNT/E12.5_CORTEX_INPUT.sorted.bam -plot $flag.png	
	done


HC


LHX2_HISTONE_HC=($(find ./ -name "*lhx2-mut*.sorted.bam" ))
for j in "${LHX2_HISTONE_HC[@]}"
	do 
	flag=$( echo "$j" | rev |cut -c 12-|rev)
	printf " plotFingerprint -b $j ./E12.5_HC_INPUT/E12.5_HC_INPUT_LHX2_MUT.sorted.bam -plot $flag.png \n"
	plotFingerprint -b $j ./E12.5_HC_INPUT/E12.5_HC_INPUT_LHX2_MUT.sorted.bam -plot $flag.png	
	done






CNT_HISTONE_HC=($(find ./ -name "*cnt*.sorted.bam" ))
for j in "${CNT_HISTONE_HC[@]}"
	do 
	printf "samtools index $j \n "
	samtools index $j 
	done

CNT_HISTONE_HC=($(find ./ -name "*cnt*.sorted.bam" ))
for j in "${CNT_HISTONE_HC[@]}"
	do 
	flag=$( echo "$j" | rev |cut -c 12-|rev)
	printf " plotFingerprint -b $j ./E12.5_HC_INPUT/E12.5_HC_INPUT_CNT.sorted.bam -plot $flag.png \n "
	plotFingerprint -b $j ./E12.5_HC_INPUT/E12.5_HC_INPUT_CNT.sorted.bam -plot $flag.png
	done


macs14 -t ./E12.5_CTX_LHX2_MUT_H3K4m1/E12.5_CTX_LHX2_MUT_H3K4m1.bam -c ./E12.5_CTX_CNT_H3K4me1/E12.5_CTX_CNT_H3K4me1.bam -f BAM -g mm -n E12.5_CORTEX_H3K4Me_CNT_LHX2_MUT_2 --single-profile --wig --space=25 --bw=300 --mfold=5,50



TAGDIRECTORY 

directory=($(ls -d e12.5*/ ))
for j in "${directory[@]}"
	do 
	printf " ./$j \n "
	cd ./$j 
	flag=($( echo "$j" |rev | cut -c 2- |rev ))
	cd ..
	printf " /home/darwin/Annu/bio_tools/bin/makeTagDirectory ./$j ./$j$flag.bam \n"
	/home/darwin/Annu/bio_tools/bin/makeTagDirectory ./$j ./$j$flag.bam
	done



/home/darwin/Annu/bio_tools/bin/makeTagDirectory ./E10.5_INPUT_REP1/ ./E10.5_INPUT_REP1/E10.5_INPUT_REP1.bam      



bw=($(find ./ -name "*.bw"))
for i in "${bw[@]}"
	do 	
	printf " cp $i /home/darwin/LHX2/BW \n"
	cp $i /home/darwin/LHX2/BW
	done



directory_ctx=($(ls -d E12.5* ))
for j in "${directory_ctx[@]}"
	do 
	cd ./$j 
	directory_ctx2=($(ls -d E12.5* ))
	for k in "${directory_ctx2[@]}"
	do 
 	flag=($( echo "$k" |rev | cut -c 1- |rev ))

	printf " /home/darwin/Annu/bio_tools/bin/makeTagDirectory ./$k ./$k/$flag.bam \n"
	/home/darwin/Annu/bio_tools/bin/makeTagDirectory ./$k ./$k/$flag.bam
	done
	cd ..
	done


findPeaks <tag directory> -style <factor | histone> -o auto -i <control tag directory>

i.e. findPeaks ERalpha-ChIP-Seq/ -style factor -o auto -i Control-ChIP-Seq/





tabel_read = read.table(txt_file_name,"\t", header=T,row.names=1)
  
  assign(SI_BUD_COMMON_MATRIC_NAME[[i]], tabel_read)
  
  converted_matrix = data.matrix(tabel_read,1:ncol(tabel_read))
  #converted_matrix_log = log(converted_matrix,1:ncol(converted_matrix))
  
  my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
  heatmap.2(converted_matrix, Colv = NA, Rowv = NA , trace = "none",col = my_palette,
            key = TRUE, key.xlab = "Log(COUNTS)", key.ylab = NULL , main = SampleB[[i]],
            xlab = "STAGES", ylab = "Log(COUNTS)")
}























