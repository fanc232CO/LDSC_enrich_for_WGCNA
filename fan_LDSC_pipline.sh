#!/bin/bash
#must in ldsc envirment.
fn_in=$1
fn_in=$(readlink -f $fn_in)
dir_in=$(echo $fn_in|awk -F'/' '{$NF=na;print}' OFS='/')
cd $dir_in
#input two columns, the first ensembl id, the second module name, no header. No grey row.
#separate by ","
dir_db=/home/fancong/database_fan
dir_ldsc=/home/fancong/softwares/ldsc

#make annotation
mkdir gene_annot
cd gene_annot
#/home/fancong/anaconda3/envs/SCAVENGE/bin/Rscript $dir_db/LDSC_enrich/scripts/biomart_loc.R $fn_in #generate gene_loc.map separated by ","
fn_loc_map=$dir_db/human_gene_loc.bed
awk 'NR==FNR{a[$1]=$2","$3","$4 }{print $1","a[$1]}' FS='\t' $fn_loc_map FS=',' $fn_in|awk -F',' 'NF==4{print}' > gene_loc.map
wait
function get_bed(){
        col=$1
        awk -F',' 'NR==FNR{a[$1]=$2"\t"$3"\t"$4;next}{if($NF==i) print a[$1]}' i=$col ./gene_loc.map $fn_in > ${col}.bed
}
for cl in $(awk -F',' '{print $2}' $fn_in|sort|uniq);do #cl--module color
        get_bed $cl
done
get_annot1(){
        col=$1 #module
        n2=$2 #chr number
        gene_bed=./${col}.bed
        bim_file=$dir_db/1000G/1000G_EUR_Phase3_plink/1000G.EUR.QC.${n2}.bim
        /home/fancong/anaconda3/envs/ldsc/bin/python2.7 $dir_ldsc/make_annot.py \
                --bed-file $gene_bed \
                --windowsize 10000 \
                --bimfile $bim_file \
                --annot-file ${n1}_${n2}.annot.gz
}
get_annot2(){
        n1=$1
        for n2 in $(seq 1 22);do
                get_annot1 $n1 $n2
        done
}
for cl in $(awk -F',' '{print $2}' $fn_in|sort|uniq);do #cl--module color
        get_annot2 $cl
done
for i in *.annot.gz;do
        gunzip $i
done
mkdir annots
mv *.annot annots/
/home/fancong/anaconda3/bin/python3.8 $dir_db/LDSC_enrich/scripts/get_chr_annot.py $fn_in
for i in $(seq 1 22);do
        gzip chr${i}.annot
done
cd ..

#calculate ldsc
mkdir ldsc
cd ldsc
fn_sumstats=$dir_db/LDSC_enrich/sumstats/lancet.sumstats.gz #change!
bn_out=$(basename $fn_sumstats .sumstats.gz)
function get_ldsc(){
        chr_n=$1
        /home/fancong/anaconda3/envs/ldsc/bin/python2.7 $dir_ldsc/ldsc.py \
        --l2 \
        --bfile $dir_db/1000G/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr_n} --ld-wind-cm 1 \
        --annot ../gene_annot/chr${chr_n}.annot.gz \
        --thin-annot \
        --out chr${chr_n} \
        #--print-snps /mnt/fanc/1000G/w_hm3.snplist
}
for chr_n in $(seq 1 22);do
        {
        echo $chr_n
        get_ldsc $chr_n
        }&
done
wait
cp ../gene_annot/chr*.annot.gz .
/home/fancong/anaconda3/envs/ldsc/bin/python2.7 $dir_ldsc/ldsc.py \
        --h2 $fn_sumstats \
        --intercept-h2 1 \
        --ref-ld-chr chr@ \
        --w-ld-chr $dir_db/1000G/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.\
        --overlap-annot \
        --frqfile-chr $dir_db/1000G/1000G_Phase3_frq/1000G.EUR.QC.\
        --out $bn_out
rm chr*.annot.gz
