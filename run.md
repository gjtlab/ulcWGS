# run

## CreateSomaticPanelOfNormals

```sh
for i in 0data_*;
do 
{ 
    out_dir=5_vcfs_0;
    ref_fa=/data/hg38/hg38.fa.gz;
    GATK=/data/apps/gatk-4.1.4.1/gatk;
    while read line; 
    do 
    { 
        sample_name=`basename $line | sed 's/.bqsr.reheader.bam//g'`; 
        mkdir -p ${out_dir}/tmps; 
        if [ ! -f ${out_dir}/${sample_name}.vcf.gz.tbi ]; 
        then 
        $GATK Mutect2 \
            --tmp-dir ${out_dir}/tmps \
            -R $ref_fa \
            -I $line \
            --max-mnp-distance 0 \
            -O ${out_dir}/${sample_name}.vcf.gz;
        fi
    }
    done < $i
} &
done

ref_fa=/data/hg38/hg38.fa.gz
GATK=/data/apps/gatk-4.1.4.1/gatk
interval=/data/GATK_reaource/hg38_v0/wgs_calling_regions.hg38.interval_list

$GATK GenomicsDBImport \
    -R $ref_fa \
    -L ${interval} \
    --merge-input-intervals TRUE \
    --genomicsdb-workspace-path pon_db \
    --reader-threads 30 \
    --batch-size 50 \
    --sample-name-map healthy_vcf.txt

$GATK CreateSomaticPanelOfNormals \
   -R /data/hg38/hg38.fa.gz \
   --germline-resource /data/GATK_reaource/somatic-hg38/af-only-gnomad.hg38.vcf.gz \
   -V gendb://nipt_pon_db \
   -O nipt_pon.vcf.gz

# add PRJNA954988_normal_tissue
$GATK GenomicsDBImport \
    -L ${interval} \
    --merge-input-intervals TRUE \
    --genomicsdb-update-workspace-path nipt_pon_db \
    --reader-threads 15 \
    --batch-size 50 \
    --sample-name-map PRJNA954988_normal_tissue_vcf.txt

mv nipt_pon_db nipt_PRJNA954988_normal_pon_db

$GATK CreateSomaticPanelOfNormals \
   -R /data/hg38/hg38.fa.gz \
   --germline-resource /data/GATK_reaource/somatic-hg38/af-only-gnomad.hg38.vcf.gz \
   -V gendb://nipt_PRJNA954988_normal_pon_db \
   -O nipt_PRJNA954988_normal_pon.vcf.gz

# add 1000g_pon
zgrep -v "#" nipt_PRJNA954988_normal_pon.vcf.gz 1000g_pon.hg38.vcf.gz | cut -d ":" -f 2- | cut -f 1-5 | awk '{print $0"\t.\t.\t."}'|sort -u -k1V -k2n > 1000G_nipt_PRJNA954988_normal_pon.txt
zgrep "#" /data/GATK_reaource/somatic-hg38/1000g_pon.hg38.vcf.gz > tmp.vcf
grep chr 1000G_nipt_PRJNA954988_normal_pon.txt >> tmp.vcf
grep -v chr 1000G_nipt_PRJNA954988_normal_pon.txt >> tmp.vcf
mv tmp.vcf 1000G_nipt_PRJNA954988_normal_pon.hg38.vcf
bgzip 1000G_nipt_PRJNA954988_normal_pon.hg38.vcf
tabix 1000G_nipt_PRJNA954988_normal_pon.hg38.vcf.gz
      1000G_nipt_PRJNA954988_normal_pon.hg38.vcf
```

## somatic variant calling based on pon

```sh
out_dir=vcfs
ref_fa=/data/hg38/hg38.fa.gz
GATK=/data/apps/gatk-4.5.0.0/gatk
pon=1000G_nipt_PRJNA954988_normal_pon.hg38.vcf.gz
snpdb=/data/GATK_reaource/somatic-hg38/af-only-gnomad.hg38.vcf.gz
while read line;
do
{
caseBam=`echo $line|awk '{print $1}'`
sample_name=`basename $caseBam | sed 's/.bqsr.bam//g'`
mkdir -p ${out_dir}/tmps
$GATK Mutect2 \
    --tmp-dir ${out_dir}/tmps \
    -R $ref_fa \
    -I $caseBam \
    --germline-resource  $snpdb \
    --panel-of-normals $pon \
    -O ${out_dir}/${sample_name}.somatic.vcf.gz
$GATK FilterMutectCalls \
    -R $ref_fa \
    -V ${out_dir}/${sample_name}.somatic.vcf.gz \
    -O ${out_dir}/${sample_name}.somatic.filtered.vcf.gz
$GATK LeftAlignAndTrimVariants \
   -R /data/hg38/hg38.fa.gz \
   -V ${out_dir}/${sample_name}.somatic.filtered.vcf.gz \
   -O ${out_dir}/${sample_name}.somatic.filtered.split.norm.vcf.gz \
   --split-multi-allelics
}
done < samples.txt
```

## merge vcf

```R
library(data.table)
library(parallel)

# data2: /data/NIPT_clean/cohort1_somatic
files = list.files("vcfs", pattern="*.bqsr.reheader.bam.somatic.filtered.split.norm.vcf.gz$", full.names=T, recursive=T)

dat = mclapply(files, function(x){
    print(x)
    d = fread(x, header=T,sep="\t",stringsAsFactors=F)
    index = d$FILTER == "PASS"
    if(sum(index)>0){
        d = d[index, ]
        if(ncol(d)==11){
            d[,10]=d[,11]
            d[,11]=NULL
        }
        colnames(d)[10] = "tumor"
        d$sample = gsub(".bqsr.reheader.bam.somatic.filtered.split.norm.vcf.gz", "", basename(x))
        d$QUAL=NULL
        d$FILTER=NULL
        return(d)
    }
}, mc.cores = 20)

dat = do.call(rbind, dat)

get_FORMAT = function(dat, format){
    dat = dat[dat$FORMAT == format, ]
    a = str_split(dat$tumor, pattern=":", simplify=T)
    colnames(a) = strsplit(format, split=":")[[1]]
    a = as.data.frame.matrix(a, stringsAsFactors=F)
    dat = cbind(dat, a)
    dat$FORMAT = NULL
    dat$tumor  = NULL
    return(dat)
}

library(stringr)
dat = lapply(unique(dat$FORMAT), function(x){
    return(get_FORMAT(dat, x))
})

dat = rbindlist(l = dat, fill = T)

info = unique(dat$INFO)
info_str = unique(unlist(strsplit(unique(str_replace_all(info, "=([^;]*)(;|$)", "=;")), split=";")))

for(i in grep("=", info_str, invert=T, v = T)){
    info = gsub(paste0(";", i, ";"), paste0(";", i, "=1;"), info)
}

get_INFO = function(info, ind, info_str){
    m = info[ind == info_str]
    mm = str_split(m, pattern="[;=]", simplify=T)
    mm = data.table(mm[,seq(2, ncol(mm), 2)])
    colnames(mm)=strsplit(info_str, split=";")[[1]]
    return(cbind(m, mm))
}

ind = str_replace_all(info, "=([^;]*)(;|$)", ";")

a = lapply(unique(ind), function(x){
    return(get_INFO(info, ind, x))
})
a = rbindlist(l=a, fill=T)
a$DP = NULL

for(i in grep("=", info_str, invert=T, v = T)){
    a$m = gsub(paste0(";", i, "=1;"), paste0(";", i, ";"), a$m)
}

dat = cbind(dat, a[match(dat$INFO, a$m),-1])
dat$INFO = NULL
dat$ID = NULL
dat = as.data.frame(dat, stringsAsFactors=F)

num_name = c("AF","DP", "PS", "CONTQ", "ECNT", "GERMQ", "MPOS", "NALOD", "NLOD", "POPAF", "SEQQ", "STRANDQ", "TLOD", "STR", "STRQ")

for(i in intersect(num_name, colnames(dat))){
    dat[, i] = as.numeric(dat[, i])
}

dp = str_split(dat$AD, pattern=",", simplify=T)
dat$RD = as.numeric(dp[,1])
dat$AD = as.numeric(dp[,2])
colnames(dat)[1] = "CHROM"
save(dat, file="71_cancers_vcfs.merge.RData")

vcf = dat[, 1:4]
vcf[, c("ID","QUAL","FILTER","INFO")] = "."
vcf = vcf[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")]
colnames(vcf)[1] = "#CHROM"
fwrite(vcf, file = "lite.vcf.gz", row.names = F, col.names = T, sep = "\t", quote = F)
```

## somatic mutations annotation

```sh
perl /data/apps/annovar/table_annovar.pl \
    lite.vcf.gz \
    -vcfinput /data/apps/annovar/humandb \
    -buildver hg38 \
    -out lite.vcf \
    -otherinfo -remove \
    -protocol refGene,avsnp150,gnomad_genome \
    -operation g,f,f -nastring .
cut -f 1-19,23-24,26-27 lite.vcf.hg38_multianno.txt | bgzip -@2 > hg38_multianno.txt.gz
rm lite.vcf.*
```

## removing germline variants

```R
library(data.table)

anno_pon = fread("nipt_PRJNA954988_normal_pon.hg38_multianno.txt.gz", header=T,sep='\t',stringsAsFactors=F)
somatic_index = which(anno_pon$gnomAD_genome_ALL == "." & !grepl("rs",anno_pon$avsnp150))

anno = anno_pon[somatic_index, ]
anno[, c("GeneDetail.refGene", "gnomAD_genome_ALL", "gnomAD_genome_AFR", "gnomAD_genome_AMR", "gnomAD_genome_ASJ", "gnomAD_genome_FIN", "gnomAD_genome_NFE", "gnomAD_genome_OTH", "gnomAD_genome_EAS", "avsnp150")] = NULL
anno$ID = paste(anno$Chr, anno$Start, anno$Ref, anno$Alt, sep="-")

ChinaMAP = fread("mbiobank_ChinaMAP_phase1_v1_reference_panel.sites.vcf.gz", header=T, sep="\t", stringsAsFactors=F)
colnames(ChinaMAP)[1] = "CHROM"
ChinaMAP = ChinaMAP[,c("CHROM", "POS","REF", "ALT", "INFO")]
ChinaMAP$INFO = gsub("AF=", "", ChinaMAP$INFO)
colnames(ChinaMAP)[5] = "AF_chinaMap"
ChinaMAP$ID = paste(ChinaMAP$CHROM, ChinaMAP$POS, ChinaMAP$REF, ChinaMAP$ALT, sep="-")

anno = anno[-which(anno$ID %in% ChinaMAP$ID), ]
anno[, c("gnomAD_genome_EAS", "avsnp150")] = NULL
anno$id = paste(anno$Otherinfo4, anno$Otherinfo5, anno$Otherinfo7, anno$Otherinfo8, sep="-")

nipt_PRJNA954988_normal_pon = fread("nipt_PRJNA954988_normal_pon.bed", header=F,sep="\t", stringsAsFactors=F)
colnames(nipt_PRJNA954988_normal_pon) = c("chr", "pos", "ref", "alt")
nipt_PRJNA954988_normal_pon$id = paste(nipt_PRJNA954988_normal_pon$chr, nipt_PRJNA954988_normal_pon$pos, nipt_PRJNA954988_normal_pon$ref, nipt_PRJNA954988_normal_pon$alt, sep="-")
nipt_PRJNA954988_normal_pon = nipt_PRJNA954988_normal_pon[which(nipt_PRJNA954988_normal_pon$id %in% anno$id), ]

tmp = nipt_PRJNA954988_normal_pon[grep(",", nipt_PRJNA954988_normal_pon$alt), ]
nipt_PRJNA954988_normal_pon = nipt_PRJNA954988_normal_pon[-grep(",", nipt_PRJNA954988_normal_pon$alt), ]

tmp_alt = str_split(tmp$alt, pattern=",")
tmp = do.call(rbind, lapply(1:nrow(tmp), function(x){
    print(x)
    return(cbind(tmp[x, -4], tmp_alt[[x]]))
}))
tmp = tmp[which(tmp$V2!="*"), ]

colnames(tmp) = c("chr", "pos", "ref", "id",  "alt")
nipt_PRJNA954988_normal_pon = rbind(nipt_PRJNA954988_normal_pon, tmp)
setkey(nipt_PRJNA954988_normal_pon, chr, pos, ref, alt)

pon_1000g = fread("1000g_pon.bed", header=F,sep="\t",stringsAsFactors=F)
colnames(pon_1000g) = c("chr", "pos", "ref", "alt")
tmp = pon_1000g[grep(",", pon_1000g$alt), ]
pon_1000g = pon_1000g[-grep(",", pon_1000g$alt), ]

tmp_alt = str_split(tmp$alt, pattern=",")
tmp = do.call(rbind, lapply(1:nrow(tmp), function(x){
    print(x)
    return(cbind(tmp[x,1:3], tmp_alt[[x]]))
}))
tmp = tmp[which(tmp$V2!="*"), ]
colnames(tmp) = c("chr", "pos", "ref", "alt")
pon_1000g = rbind(pon_1000g, tmp)
setkey(pon_1000g, chr, pos, ref, alt)

nipt_PRJNA954988_normal_pon_rm1000g = nipt_PRJNA954988_normal_pon[!pon_1000g]

anno = anno[anno$id %in% nipt_PRJNA954988_normal_pon_rm1000g$id, ]

bases = c("A", "T", "C", "G")
names(bases) = c("T", "A", "G", "C")

### triBase
get_triBase = function(mut, genome_version = "hg38"){
    colnames(mut) = c("chr", "pos", "ref", "alt")
    mat <- mut[-which(mut$ref %in% c("A", "T", "C", "G") & mut$alt %in% c("A", "T", "C", "G")), ]
    mat$ref1 = mat$ref
    mat$alt1 = mat$alt
    mut <- mut[which(mut$ref %in% c("A", "T", "C", "G") & mut$alt %in% c("A", "T", "C", "G")), ]
    if(genome_version == "hg38"){
        hsapiens = BSgenome.Hsapiens.UCSC.hg38::Hsapiens
    }else if(genome_version == "hg19"){
        hsapiens = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
    }
    mut_positive = mut[which(mut$ref %in% c("C", "T"))]
    mut_negative = mut[which(mut$ref %in% c("G", "A"))]
    mut_positive$context = BSgenome::getSeq(hsapiens, mut_positive$chr, mut_positive$pos - 1, mut_positive$pos + 1, as.character = T)
    mut_negative$context = BSgenome::getSeq(hsapiens, mut_negative$chr, mut_negative$pos - 1, mut_negative$pos + 1, as.character = T, strand="-")
    mut_positive = mut_positive[which(mut_positive$ref == substr(mut_positive$context, 2, 2)), ]
    mut_positive$ref1 = mut_positive$ref
    mut_positive$alt1 = mut_positive$alt
    mut_negative$ref1 = bases[mut_negative$ref]
    mut_negative$alt1 = bases[mut_negative$alt]
    mut_negative = mut_negative[which(mut_negative$ref1 == substr(mut_negative$context, 2, 2)), ]
    mut = rbind(mut_positive, mut_negative)
    mut$var = paste0(mut$ref1, ">", mut$alt1)
    mut$tricontext = paste(substr(mut$context, 1, 1), "[", mut$var, "]", substr(mut$context, 3, 3), sep = "")
    mut <- mut[, c("chr", "pos", "ref", "alt", "ref1", "alt1", "tricontext")]T
    mat$tricontext <- ""
    mut <- rbind(mut, mat)
    return(mut)
}
mut = get_triBase(anno[,c(1:2,4:5)])
mut$ID = paste(mut$chr, mut$pos, mut$ref, mut$alt, sep="-")
anno$tricontext = mut$tricontext[match(anno$ID, mut$ID)]

fwrite(anno, file="nipt_PRJNA954988_normal_pon.hg38_multianno.rmSNP150.rmChinaMap.rm1000gPon.txt.gz", row.names=F,col.names=T, sep="\t",quote=F)

###################

anno_raw = fread("vcfs_merge.hg38_multianno.txt.gz", header=T,sep='\t',stringsAsFactors=F)
somatic_index = which(anno_raw$gnomAD_genome_ALL == "." & !grepl("rs", anno_raw$avsnp150))
anno_raw[, c("GeneDetail.refGene", "gnomAD_genome_ALL", "gnomAD_genome_AFR", "gnomAD_genome_AMR", "gnomAD_genome_ASJ", "gnomAD_genome_FIN", "gnomAD_genome_NFE", "gnomAD_genome_OTH", "gnomAD_genome_EAS", "avsnp150")] = NULL
anno_raw = anno_raw[somatic_index, ]
setkey(anno_raw, Chr, Start, End, Ref, Alt)

anno = fread("nipt_PRJNA954988_normal_pon.hg38_multianno.rmSNP150.rmChinaMap.rm1000gPon.txt.gz", header=T,sep="\t",stringsAsFactors=F)  
setkey(anno, Chr, Start, End, Ref, Alt)

a = anno_raw[anno, nomatch = 0]
a[, c("i.Func.refGene", "i.Gene.refGene", "i.ExonicFunc.refGene", "i.AAChange.refGene", "i.Otherinfo4", "i.Otherinfo5", "i.Otherinfo7", "i.Otherinfo8", "id", "ID")] = NULL

fwrite(a, file="nipt_pon.hg38_multianno.rmSNP150.rmChinaMap.rm1000gPon.txt.gz", row.names=F,col.names=T, sep="\t",quote=F)

load("vcfs.merge.RData")
dat$ID = paste(dat$CHROM, dat$POS, dat$REF, dat$ALT, sep="-")
dat = dat[which(dat$ID %in% anno$id), ]

pm = dat[dat$REF %in% c("A", "T", "G", "C") & dat$ALT %in% c("A", "T", "G", "C"), c("sample", "CHROM", "POS", "REF", "ALT")]
fwrite(pm, file="vcfs.merge.pm.gz", row.names=F,col.names=F,sep="\t",quote=F)
```

## somatic signatures

```sh
mkdir vcfs_somatic_deconstructSigs_out
cd vcfs_somatic_deconstructSigs_out
Rscript deconstructSigs/deconstructSigs.R -i vcfs.merge.pm.gz -o sig_output -p 20 -r COSMIC_v3.3.1_SBS_GRCh38.refSignature.txt -c 0 -b hg38
```

```R
load("sig_output.RData")
sigMat = as.data.frame(t(sigs.weights))
write.table(sigMat, file="sigMat.tsv", r=T,c=T,sep="\t",quote=F)
```
