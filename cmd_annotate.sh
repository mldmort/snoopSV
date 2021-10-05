#!/bin/bash

CHR="$1"
SAMPLE_DEPTH_DIR="$2"
SAMPLE_BAM_FILE="$3"
VCF_IN="$4"
ANNOT_OUT="$5"
TMP_DIR="$6"

echo -e "chrom\tpos\tend\tid\tsvtype\tchr2\tend2\trv\tdp_samples\tdp_fr_samples\tmapq_all_samples\tsamples" > $ANNOT_OUT

declare -A sample_dep_dict
declare -A sample_dep_dir_dict
declare -A sample_gen_dict
declare -A sample_bam_dict
declare -A main_chr_dict
declare -a sample_list


while read -r line; do
	sample=$(echo "$line" | awk '{print $1}')
	dir=$(echo "$line" | awk '{print $2}')
	sample_dep_dict[$sample]=$(ls $dir/*$CHR.bed)
	sample_gen_dict[$sample]=$(ls $dir/*$CHR.genome)
	sample_dep_dir_dict[$sample]=$dir
done < $SAMPLE_DEPTH_DIR

while read -r line; do
	sample=$(echo "$line" | awk '{print $1}')
	bam=$(echo "$line" | awk '{print $2}')
	sample_bam_dict[$sample]=$bam
done < $SAMPLE_BAM_FILE

for sample in $(bcftools query -l $VCF_IN); do
	sample_list+=( $sample )
done

for i in {1..24}; do
	CHR=chr$i
	if [ $i == 23 ]; then
		CHR=chrX
	elif [ $i == 24 ]; then
		CHR=chrY
	fi
	main_chr_dict[$CHR]=1
done

#echo "keys:" ${!sample_dep_dict[@]} " vals: " ${sample_dep_dict[@]}
#echo "keys:" ${!sample_gen_dict[@]} " vals: " ${sample_gen_dict[@]}
#echo "keys:" ${!sample_bam_dict[@]} " vals: " ${sample_bam_dict[@]}
#echo "keys:" ${!sample_list[@]} " vals: " ${sample_list[@]}

l_thr=100000
while read -r line; do
	#echo "$line"

	svtype=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print $5}')

	samples=""
	sep=""
	for sample in ${sample_list[@]}; do
		samples=$samples$sep$sample
		sep=","
	done

	if [ $svtype != "TRA" ]; then
		this_line=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print $1, $2, $3}')
		chr_1=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print $1}')
		chr_2=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print $6}')
		
		pos=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print $2}')
		end=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print $3}')
		svlen=$(expr $end - $pos)
		fr_len=$(expr $svlen / 3 )
		if [ $fr_len -gt $l_thr ]; then
		    fr_len=$l_thr
		fi  
		pos_fr_1=$(expr $pos - $fr_len)
		pos_fr_2=$pos
		end_fr_1=$end
		end_fr_2=$(expr $end + $fr_len)
		if [ $pos_fr_1 -lt 0 ]; then
		    pos_fr_1=0
		fi  
		region_fr=$chr_1:$pos_fr_1-$pos_fr_2" "$chr_1:$end_fr_1-$end_fr_2
		#echo "region_fr: $region_fr"
		
		#echo "svtype: $svtype"
		#echo "pos_fr_1: $pos_fr_1"
		#echo "pos_fr_2: $pos_fr_2"
		#echo "end_fr_1: $end_fr_1"
		#echo "end_fr_2: $end_fr_2"


		mapq_samples_all=""
		sep=""
		for sample in ${sample_list[@]}; do
			file_bam=${sample_bam_dict[$sample]}
			this_mapq=$(samtools view $file_bam $region_fr | sort -u -T $TMP_DIR | \
						awk 'BEGIN{FS="\t";OFS="\t"; q_sum_all=0; n_all=0; } \
							{q_sum_all+=$5; n_all+=1; } \
							END{if(n_all>0){q_all=q_sum_all/n_all}else{q_all=0}; \
								print q_all}')
			mapq_samples_all=$mapq_samples_all$sep$this_mapq
			sep=","
		done
		#echo "mapq all: $mapq_samples_all"

		depth_samples=""
		sep=""
		for sample in ${sample_list[@]}; do
			file_depth=${sample_dep_dict[$sample]}
			file_genome=${sample_gen_dict[$sample]}
			this_depth=$(bedtools intersect -sorted -g $file_genome -a <(echo "$this_line") -b $file_depth -wo | \
			awk 'BEGIN{FS="\t";OFS="\t";sum_dep=0;sum_ovr=0} \
			    {sum_dep+=$7 * $8; sum_ovr+=$8;} \
			    END{if(sum_ovr>0){print sum_dep/sum_ovr}else{print 0}}')
			depth_samples=$depth_samples$sep$this_depth
			sep=","
		done
		#echo "depth_samples: $depth_samples"

		line_fr_1=$(echo -e "$chr_1\t$pos_fr_1\t$pos_fr_2")
		line_fr_2=$(echo -e "$chr_1\t$end_fr_1\t$end_fr_2")
		#echo "line_fr_1: $line_fr_1"
		#echo "line_fr_2: $line_fr_2"
		 
		depth_flank_samples=""
		sep=""
		for sample in ${sample_list[@]}; do
			file_depth=${sample_dep_dict[$sample]}
			file_genome=${sample_gen_dict[$sample]}
			this_depth_fr_1=$(bedtools intersect -sorted -g $file_genome -a <(echo "$line_fr_1") -b $file_depth -wo | \
			awk 'BEGIN{FS="\t";OFS="\t";sum_dep=0;sum_ovr=0} \
				{sum_dep+=$7 * $8; sum_ovr+=$8;} \
				END{if(sum_ovr>0){print sum_dep/sum_ovr}else{print 0}}')
			this_depth_fr_2=$(bedtools intersect -sorted -g $file_genome -a <(echo "$line_fr_2") -b $file_depth -wo | \
			awk 'BEGIN{FS="\t";OFS="\t";sum_dep=0;sum_ovr=0} \
				{sum_dep+=$7 * $8; sum_ovr+=$8;} \
				END{if(sum_ovr>0){print sum_dep/sum_ovr}else{print 0}}')
			depth_flank_samples=$depth_flank_samples$sep$this_depth_fr_1:$this_depth_fr_2
			sep=","
		done
		#echo "depth_flank_samples: $depth_flank_samples"
	elif [ $svtype == "TRA" ]; then
		fr_len=5000
		#'%CHROM\t%POS0\t%END\t%ID\t%SVTYPE\t%CHR2\t%AVG_END\n'

		chr_1=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print $1}')
		pos_1=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print $2}')
		end_1=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print $2+1}')
		chr_2=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print $6}')
		pos_2=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print int($7)}')
		end_2=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print int($7)+1}')

		this_line_1=$(echo -e "$chr_1\t$pos_1\t$end_1")
		this_line_2=$(echo -e "$chr_2\t$pos_2\t$end_2")
		#echo "this_line_1: $this_line_1"
		#echo "this_line_2: $this_line_2"

		pos_fr_1=$(expr $pos_1 - $fr_len)
		pos_fr_2=$(expr $pos_1 + $fr_len)
		end_fr_1=$(expr $pos_2 - $fr_len)
		end_fr_2=$(expr $pos_2 + $fr_len)
		if [ $pos_fr_1 -lt 0 ]; then
			pos_fr_1=0
		fi
		if [ $end_fr_1 -lt 0 ]; then
			end_fr_1=0
		fi
		region_fr=$chr_1:$pos_fr_1-$pos_fr_2" "$chr_2:$end_fr_1-$end_fr_2
		#echo "region_fr: $region_fr"

		mapq_samples_all=""
		sep=""
		for sample in ${sample_list[@]}; do
			if [ ${main_chr_dict[$chr_2]} ]; then
				file_bam=${sample_bam_dict[$sample]}
				this_mapq=$(samtools view $file_bam $region_fr | sort -u -T $TMP_DIR | \
							awk 'BEGIN{FS="\t";OFS="\t"; q_sum_all=0; n_all=0; } \
								{q_sum_all+=$5; n_all+=1; } \
								END{if(n_all>0){q_all=q_sum_all/n_all}else{q_all=0}; \
									print q_all}')
			else
				this_mapq="."
			fi
			mapq_samples_all=$mapq_samples_all$sep$this_mapq
			sep=","
		done
		#echo "mapq all: $mapq_samples_all"

		depth_samples=""
		sep=""
		for sample in ${sample_list[@]}; do
			if [ ${main_chr_dict[$chr_2]} ]; then
				file_depth_dir=${sample_dep_dir_dict[$sample]}
				file_depth_1=$(ls $file_depth_dir/*$chr_1.bed)
				file_genome_1=$(ls $file_depth_dir/*$chr_1.genome)
				file_depth_2=$(ls $file_depth_dir/*$chr_2.bed)
				file_genome_2=$(ls $file_depth_dir/*$chr_2.genome)
				#echo "file_depth_1: $file_depth_1"
				#echo "file_depth_2: $file_depth_2"
				#echo "file_genome_1: $file_genome_1"
				#echo "file_genome_2: $file_genome_2"

				this_depth_1=$(bedtools intersect -sorted -g $file_genome_1 -a <(echo "$this_line_1") -b $file_depth_1 -wo | awk 'BEGIN{FS="\t";OFS="\t"}{print $7}')
				this_depth_2=$(bedtools intersect -sorted -g $file_genome_2 -a <(echo "$this_line_2") -b $file_depth_2 -wo | awk 'BEGIN{FS="\t";OFS="\t"}{print $7}')
			else
				this_depth_1="."
				this_depth_2="."
			fi
			depth_samples=$depth_samples$sep$this_depth_1:$this_depth_2
			sep=","
		done
		#echo "depth_samples: $depth_samples"

		depth_flank_samples="."
		#echo "depth_flank_samples: $depth_flank_samples"
	fi
	
	echo -e "$line\t$depth_samples\t$depth_flank_samples\t$mapq_samples_all\t$samples" >> $ANNOT_OUT

done < <(bcftools query -e 'SKIP_REGION==1 || SKIP_TR==1' -f '%CHROM\t%POS0\t%END\t%ID\t%SVTYPE\t%CHR2\t%AVG_END\t[%RV,]\n' $VCF_IN | awk 'BEGIN{FS="\t";OFS="\t"}{if($5=="INS"){$2-=50;$3+=49}; $NF=substr($NF, 1, length($NF)-1); print $0}')

