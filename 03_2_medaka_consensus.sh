regions=`cut -f 1 assembly_info.txt | sed -n '2,3000p'`

medaka consensus calls_to_draft.bam.bam ./split_result/sub_1.hdf --model r1041_e82_400bps_hac_v4.2.0 --batch 100 --threads 2 --region ${regions}