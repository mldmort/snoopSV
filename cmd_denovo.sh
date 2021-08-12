#!/bin/bash

BED_IN="$1"
BED_OUT="$2"
SQ_THR="$3"

REF=/expanse/projects/sebat1/genomicsdataanalysis/00001/00001_data/00001_reference-genome/Homo_sapiens_assembly38.fasta
DIR_ILL=/expanse/projects/sebat1/genomicsdataanalysis/REACH_MSSNG/request_202001

#echo "BED_IN: $BED_IN"
#echo "SQ_THR: $SQ_THR"
#echo "BED_OUT: $BED_OUT"

col_sq=$(head -n 1 $BED_IN | awk 'BEGIN{FS="\t";OFS="\t"}{for(i=1;i<=NF;i++){if($i=="SQ"){print i}}}')
col_dn=$(head -n 1 $BED_IN | awk 'BEGIN{FS="\t";OFS="\t"}{for(i=1;i<=NF;i++){if($i=="de-novo"){print i}}}')
col_n_sam_rv=$(head -n 1 $BED_IN | awk 'BEGIN{FS="\t";OFS="\t"}{for(i=1;i<=NF;i++){if($i=="N_sample_rv"){print i}}}')

col_chrom=$(head -n 1 $BED_IN | awk 'BEGIN{FS="\t";OFS="\t"}{for(i=1;i<=NF;i++){if($i=="chrom"){print i}}}')
col_pos=$(head -n 1 $BED_IN | awk 'BEGIN{FS="\t";OFS="\t"}{for(i=1;i<=NF;i++){if($i=="pos"){print i}}}')
col_end=$(head -n 1 $BED_IN | awk 'BEGIN{FS="\t";OFS="\t"}{for(i=1;i<=NF;i++){if($i=="end"){print i}}}')
col_svtype=$(head -n 1 $BED_IN | awk 'BEGIN{FS="\t";OFS="\t"}{for(i=1;i<=NF;i++){if($i=="svtype"){print i}}}')
col_sample=$(head -n 1 $BED_IN | awk 'BEGIN{FS="\t";OFS="\t"}{for(i=1;i<=NF;i++){if($i=="sample"){print i}}}')
col_mom=$(head -n 1 $BED_IN | awk 'BEGIN{FS="\t";OFS="\t"}{for(i=1;i<=NF;i++){if($i=="mom"){print i}}}')
col_dad=$(head -n 1 $BED_IN | awk 'BEGIN{FS="\t";OFS="\t"}{for(i=1;i<=NF;i++){if($i=="dad"){print i}}}')

echo -e "chrom\tpos\tend\tsvtype\tsample\tmom\tdad\tdepth_sam\tdepth_mom\tdepth_dad\tdepth_sam_fr\tdepth_mom_fr\tdepth_dad_fr\tN_PAIR_SAM\tN_PAIR_MOM\tN_PAIR_DAD\tjumps_sam\tjumps_mom\tjumps_dad" > $BED_OUT

l_thr=1000

while read -r line; do

	this_line=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print $1, $2, $3}')
	chrom=$(echo "$line" | awk -v col_chrom="$col_chrom" 'BEGIN{FS="\t";OFS="\t"}{print $col_chrom}')
	pos=$(echo "$line" | awk -v col_pos="$col_pos" 'BEGIN{FS="\t";OFS="\t"}{print $col_pos}')
	end=$(echo "$line" | awk -v col_end="$col_end" 'BEGIN{FS="\t";OFS="\t"}{print $col_end}')
	svtype=$(echo "$line" | awk -v col_svtype="$col_svtype" 'BEGIN{FS="\t";OFS="\t"}{print $col_svtype}')
	sample=$(echo "$line" | awk -v col_sample="$col_sample" 'BEGIN{FS="\t";OFS="\t"}{print $col_sample}')
	mom=$(echo "$line" | awk -v col_mom="$col_mom" 'BEGIN{FS="\t";OFS="\t"}{print $col_mom}')
	dad=$(echo "$line" | awk -v col_dad="$col_dad" 'BEGIN{FS="\t";OFS="\t"}{print $col_dad}')
	SAMPLE=$(echo $sample | awk '{print toupper($1)}')
	MOM=$(echo $mom | awk '{print toupper($1)}')
	DAD=$(echo $dad | awk '{print toupper($1)}')

	fl_reg=$(echo "$line" | \
		awk -v l_thr="$l_thr" \
			-v col_svtype="$col_svtype" \
			'BEGIN{FS="\t";OFS="\t"}\
			{\
				if($col_svtype == "INS"){\
					print $2-100, $2, $3, $3+100;\
				}\
				else{\
					svlen = $3 - $2;\
					l = int(svlen / 3);\
					if (l > l_thr) {l = l_thr};\
					print $2-l, $2, $3, $3+l;\
				}\
			}')

	pos_fr_1=$(echo "$fl_reg" | awk 'BEGIN{FS="\t";OFS="\t"}{print $1}')
	pos_fr_2=$(echo "$fl_reg" | awk 'BEGIN{FS="\t";OFS="\t"}{print $2}')
	end_fr_1=$(echo "$fl_reg" | awk 'BEGIN{FS="\t";OFS="\t"}{print $3}')
	end_fr_2=$(echo "$fl_reg" | awk 'BEGIN{FS="\t";OFS="\t"}{print $4}')

	REG=$(echo "$chrom:$pos-$end")
	FL_REG_1=$(echo "$chrom:$pos_fr_1-$pos_fr_2")
	FL_REG_2=$(echo "$chrom:$end_fr_1-$end_fr_2")
	#echo "REG: $REG"
	#echo "FL_REG_1: $FL_REG_1"
	#echo "FL_REG_2: $FL_REG_2"

	cram_sam=$(ls $DIR_ILL/$SAMPLE\_recal.cram)
	cram_mom=$(ls $DIR_ILL/$MOM\_recal.cram)
	cram_dad=$(ls $DIR_ILL/$DAD\_recal.cram)

	depth_sam=$(samtools depth -r $REG --reference $REF $cram_sam | awk 'BEGIN{FS="\t";OFS="\t";sum_dep=0}{sum_dep+=$3}END{if(NR>0){print sum_dep/NR}else{print 0}}')
	depth_mom=$(samtools depth -r $REG --reference $REF $cram_mom | awk 'BEGIN{FS="\t";OFS="\t";sum_dep=0}{sum_dep+=$3}END{if(NR>0){print sum_dep/NR}else{print 0}}')
	depth_dad=$(samtools depth -r $REG --reference $REF $cram_dad | awk 'BEGIN{FS="\t";OFS="\t";sum_dep=0}{sum_dep+=$3}END{if(NR>0){print sum_dep/NR}else{print 0}}')

	depth_sam_fr_1=$(samtools depth -r $FL_REG_1 --reference $REF $cram_sam | awk 'BEGIN{FS="\t";OFS="\t";sum_dep=0}{sum_dep+=$3}END{if(NR>0){print sum_dep/NR}else{print 0}}')
	depth_mom_fr_1=$(samtools depth -r $FL_REG_1 --reference $REF $cram_mom | awk 'BEGIN{FS="\t";OFS="\t";sum_dep=0}{sum_dep+=$3}END{if(NR>0){print sum_dep/NR}else{print 0}}')
	depth_dad_fr_1=$(samtools depth -r $FL_REG_1 --reference $REF $cram_dad | awk 'BEGIN{FS="\t";OFS="\t";sum_dep=0}{sum_dep+=$3}END{if(NR>0){print sum_dep/NR}else{print 0}}')

	depth_sam_fr_2=$(samtools depth -r $FL_REG_2 --reference $REF $cram_sam | awk 'BEGIN{FS="\t";OFS="\t";sum_dep=0}{sum_dep+=$3}END{if(NR>0){print sum_dep/NR}else{print 0}}')
	depth_mom_fr_2=$(samtools depth -r $FL_REG_2 --reference $REF $cram_mom | awk 'BEGIN{FS="\t";OFS="\t";sum_dep=0}{sum_dep+=$3}END{if(NR>0){print sum_dep/NR}else{print 0}}')
	depth_dad_fr_2=$(samtools depth -r $FL_REG_2 --reference $REF $cram_dad | awk 'BEGIN{FS="\t";OFS="\t";sum_dep=0}{sum_dep+=$3}END{if(NR>0){print sum_dep/NR}else{print 0}}')

	if [ $svtype == "INS" ]; then
		REG_INS=$(echo "$chrom:$pos_fr_1-$end_fr_2")
		#echo "REG_INS: $REG_INS"
		N_PAIR_SAM=$(samtools view -F 12 --reference $REF $cram_sam $REG_INS | awk 'BEGIN{FS="\t";OFS="\t"}$7!="="{print $7}' | wc -l)
		N_PAIR_MOM=$(samtools view -F 12 --reference $REF $cram_mom $REG_INS | awk 'BEGIN{FS="\t";OFS="\t"}$7!="="{print $7}' | wc -l)
		N_PAIR_DAD=$(samtools view -F 12 --reference $REF $cram_dad $REG_INS | awk 'BEGIN{FS="\t";OFS="\t"}$7!="="{print $7}' | wc -l)

		jumps_sam=$(samtools depth -r $REG_INS --reference $REF $cram_sam | \
		awk 'BEGIN{FS="\t";OFS="\t";mean=0;std=0;N=0;split("",diff_dic,"")}\
			{\
				this_cov = $3;\
				if(NR>1){\
					diff = this_cov - prev_cov;\
				}\
				else{\
					diff = 0;\
				}\
				prev_cov = $3;\
				mean+=diff; std+=diff*diff; N+=1;\
				diff_dic[$2] = diff;\
			}\
			END{\
				mean = mean/N;\
				std = std/N - mean*mean;\
				jumps=""; sep="";\
				for (key in diff_dic) {\
					this_diff = sqrt(diff_dic[key]*diff_dic[key]);\
					if (this_diff > std*3) {\
						jumps=jumps""sep""key"|"this_diff;\
						sep=",";\
					}\
				}\
				if (jumps == "") {jumps = "."};\
				print jumps;\
			}')
		jumps_mom=$(samtools depth -r $REG_INS --reference $REF $cram_mom | \
		awk 'BEGIN{FS="\t";OFS="\t";mean=0;std=0;N=0;split("",diff_dic,"")}\
			{\
				this_cov = $3;\
				if(NR>1){\
					diff = this_cov - prev_cov;\
				}\
				else{\
					diff = 0;\
				}\
				prev_cov = $3;\
				mean+=diff; std+=diff*diff; N+=1;\
				diff_dic[$2] = diff;\
			}\
			END{\
				mean = mean/N;\
				std = std/N - mean*mean;\
				jumps=""; sep="";\
				for (key in diff_dic) {\
					this_diff = sqrt(diff_dic[key]*diff_dic[key]);\
					if (this_diff > std*3) {\
						jumps=jumps""sep""key"|"this_diff;\
						sep=",";\
					}\
				}\
				if (jumps == "") {jumps = "."};\
				print jumps;\
			}')
		jumps_dad=$(samtools depth -r $REG_INS --reference $REF $cram_dad | \
		awk 'BEGIN{FS="\t";OFS="\t";mean=0;std=0;N=0;split("",diff_dic,"")}\
			{\
				this_cov = $3;\
				if(NR>1){\
					diff = this_cov - prev_cov;\
				}\
				else{\
					diff = 0;\
				}\
				prev_cov = $3;\
				mean+=diff; std+=diff*diff; N+=1;\
				diff_dic[$2] = diff;\
			}\
			END{\
				mean = mean/N;\
				std = std/N - mean*mean;\
				jumps=""; sep="";\
				for (key in diff_dic) {\
					this_diff = sqrt(diff_dic[key]*diff_dic[key]);\
					if (this_diff > std*3) {\
						jumps=jumps""sep""key"|"this_diff;\
						sep=",";\
					}\
				}\
				if (jumps == "") {jumps = "."};\
				print jumps;\
			}')
		#echo "jumps_sam: $jumps_sam"
		#echo "jumps_mom: $jumps_mom"
		#echo "jumps_dad: $jumps_dad"
	else
		N_PAIR_SAM="."
		N_PAIR_MOM="."
		N_PAIR_DAD="."
		jumps_sam="."
		jumps_mom="."
		jumps_dad="."
	fi

	echo -e "$chrom\t$pos\t$end\t$svtype\t$sample\t$mom\t$dad\t$depth_sam\t$depth_mom\t$depth_dad\t$depth_sam_fr_1|$depth_sam_fr_2\t$depth_mom_fr_1|$depth_mom_fr_2\t$depth_dad_fr_1|$depth_dad_fr_2\t$N_PAIR_SAM\t$N_PAIR_MOM\t$N_PAIR_DAD\t$jumps_sam\t$jumps_mom\t$jumps_dad" >> $BED_OUT

done < <(\
awk -v col_sq="$col_sq" \
	-v col_dn="$col_dn" \
	-v col_n_sam_rv="$col_n_sam_rv" \
	-v sq_thr="$SQ_THR" \
	'BEGIN{FS="\t";OFS="\t"} $1!="chrom" && $col_dn==1 && $col_n_sam_rv==1 && $col_sq>=sq_thr' $BED_IN)
