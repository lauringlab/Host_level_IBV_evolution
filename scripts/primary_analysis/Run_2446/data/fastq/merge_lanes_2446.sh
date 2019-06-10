ls -1 *R1*.gz | awk -F '_' '{print $1 "_" $2 "_" $3}' | sort | uniq > IDs.txt

for i in `cat ./IDs.txt`; do 
	var=$(echo $i | awk -F '_' '{print $1 "_" $3}')
	echo $var
	cat "$i"_L001_R1_001.fastq.gz "$i"_L002_R1_001.fastq.gz "$i"_L003_R1_001.fastq.gz "$i"_L004_R1_001.fastq.gz > "$var"_L001_R1_001.fastq.gz
done

for i in `cat ./IDs.txt`; do
        var=$(echo $i | awk -F '_' '{print $1 "_" $3}')
        echo $var
        cat "$i"_L001_R2_001.fastq.gz "$i"_L002_R2_001.fastq.gz "$i"_L003_R2_001.fastq.gz "$i"_L004_R2_001.fastq.gz > "$var"_L001_R2_001.fastq.gz
done

rm IDs.txt
# echo rm *-*
