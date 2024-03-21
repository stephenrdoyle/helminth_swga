# SWGA development

## Reference sequences

```bash



# reference using all STH but not Ttrichuria
cat ancylostoma_ceylanicum.PRJNA231479.WBPS16.genomic.fa ascaris_lumbricoides.PRJEB4950.WBPS16.genomic.fa ancylostoma_duodenale.PRJNA72581.WBPS16.genomic.fa strongyloides_stercoralis.PRJEB528.WBPS16.genomic.fa necator_americanus.PRJNA72135.WBPS16.genomic.fa > sth_minus_ttrichuria.fa
```



## swga.pl
- script from https://doi.org/10.1534/genetics.114.165498
- used for malaria SWGA design, ie. https://doi.org/10.1186/s12936-016-1641-7  

```bash
# working directory
/nfs/users/nfs_s/sd21/lustre118_link/STH/SWGA

cp ../../trichuris_trichiura/01_REF/trichuris_trichiura.fa .

wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.fna.gz


# run swga.pl targeting Trichuris with human as contaminant
#--- target: trichuris genome trichuris_trichiura.PRJEB535.WBPS16.genomic.fa
#--- offtarget: human genome GRCh38_latest_genomic.fna
#--- primer length = 12
# Tm limit = 30
#--- number of hits = 100
#--- output filename =

bsub.py 5 swga_tt_v_hs "perl swga.pl trichuris_trichiura.PRJEB535.WBPS16.genomic.fa GRCh38_latest_genomic.fna 12 30 100 ttrichiura_v_human_swga-out"

bsub.py 5 swga_tt2_v_hs "perl swga.pl trichuris_trichiura.fa GRCh38_latest_genomic.fna 12 30 100 ttrichiura2_v_human_swga-out"






cd /nfs/users/nfs_s/sd21/lustre118_link/STH/SWGA_TRICHURIS_V_STH

ln -s ../REF_GENOMES/trichuris_trichiura.fa
ln -s ../REF_GENOMES/sth_minus_ttrichuria.fa

bsub.py 5 swga_tt_v_sth "perl /nfs/users/nfs_s/sd21/lustre118_link/STH/SWGA/swga.pl trichuris_trichiura.fa sth_minus_ttrichuria.fa 12 30 100 ttrichiura_v_sth_swga-out"

cat ttrichiura_v_sth_swga-out | sort -k4,4n | tail -n20
```
- ranked 20 list of primers, hits on target, hits on contaminant, target/contaminant enrichment, Tm

|--------------|-------------|------------------|--------------------------|----|
| primer       | target_hits | contaminant_hits | target_contaminant_ratio | Tm |
|--------------|-------------|------------------|--------------------------|----|
| AAGTAAACTTCA | 894         | 448              | 1.99553571428571         | 30 |
| TTCCGTAATAAA | 635         | 317              | 2.00315457413249         | 30 |
| GCTTTAAACTTT | 600         | 261              | 2.29885057471264         | 30 |
| CAAATTGTTAGA | 558         | 240              | 2.325                    | 30 |
| TCCGTAATAAAA | 546         | 219              | 2.49315068493151         | 30 |
| TTTAATTACAGA | 741         | 290              | 2.5551724137931          | 28 |
| TAGTTTAAAGTG | 537         | 195              | 2.75384615384615         | 30 |
| TATCTGTAATTA | 765         | 274              | 2.79197080291971         | 28 |
| TTAATTACAGAT | 676         | 233              | 2.90128755364807         | 28 |
| TGTTAAAAGACA | 634         | 196              | 3.23469387755102         | 30 |
| GTTTAATTACAG | 670         | 176              | 3.80681818181818         | 30 |
| AAGTTTACTTAC | 930         | 227              | 4.09691629955947         | 30 |
| TTACAGATACAT | 854         | 205              | 4.16585365853659         | 30 |
| AGTTTAAAGTGA | 1134        | 259              | 4.37837837837838         | 30 |
| TTTACTTACAGA | 1004        | 205              | 4.89756097560976         | 30 |
| TAAGTAAACTTC | 904         | 177              | 5.10734463276836         | 30 |
| AGTTTACTTACA | 1009        | 197              | 5.12182741116751         | 30 |
| TTACTTACAGAT | 960         | 173              | 5.54913294797688         | 30 |
| ATCCGTAATAAA | 1170        | 187              | 6.25668449197861         | 30 |
| TATCTGTAAGTA | 963         | 134              | 7.1865671641791          | 30 |

- the enrichment values seem relative modest, but the contaminant is ALL other STH species whcih add up to about 1.2 Gb, compared with ~80 Mb trichuris, so the enrichment is far greater.
- will order these top 20 primers, adding a phosphorothioate bonds between the two most 3′ nucleotides to prevent primer degradation by phi29, ie  TATCTGTAAG\*T\*A

```bash
# extract the coordinates of the hits
cat ttrichiura_v_sth_swga-out | sort -k4,4n | cut -f1 | tail -n20 | while read SEQ; do
     fastaq search_for_seq trichuris_trichiura.fa trichuris_trichiura.${SEQ}.coords ${SEQ};
done

cat *.coords | sort -k1,1 -k2,2n > trichuris_trichiura.swga_primer_hits.txt
cat trichuris_trichiura.swga_primer_hits.txt | awk '{print $1,$2,$2+1,$3}' OFS="\t" >trichuris_trichiura.swga_primer_hits.bed

cat trichuris_trichiura.fa.fai | awk '{print $1,$2}' OFS="\t" > trichuris_trichiura.genome

bedtools makewindows -g trichuris_trichiura.genome -w 50000 > trichuris_trichiura.50k_window.bed

bedtools coverage -a trichuris_trichiura.50k_window.bed -b trichuris_trichiura.swga_primer_hits.bed > trichuris_trichiura.swga_primer_hits.50k.coverage



```


### SWGA - Haemonchus vs tcircumcincta
```bash
cd /nfs/users/nfs_s/sd21/lustre118_link/STH/SWGA_HAEM_V_TCIRC

wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS16.genomic.fa.gz

wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/teladorsagia_circumcincta/PRJNA72569/teladorsagia_circumcincta.PRJNA72569.WBPS16.genomic.fa.gz

gunzip *

bsub.py 5 swga_haem_v_tcirc "perl /nfs/users/nfs_s/sd21/lustre118_link/STH/SWGA/swga.pl haemonchus_contortus.PRJEB506.WBPS16.genomic.fa teladorsagia_circumcincta.PRJNA72569.WBPS16.genomic.fa 12 30 100 haem_v_tcirc_swga-out"

cat haem_v_tcirc_swga-out | sort -k4,4n | tail -n20

# extract the coordinates of the hits
cat haem_v_tcirc_swga-out | sort -k4,4n | cut -f1 | tail -n50 | while read SEQ; do
     fastaq search_for_seq haemonchus_contortus.PRJEB506.WBPS16.genomic.fa ${SEQ}.coords ${SEQ};
done

cat *.coords | sort -k1,1 -k2,2n > swga_primer_hits.txt
cat swga_primer_hits.txt | awk -F "\t" '{print $1,$2,$2+1,$3}' OFS="\t" > swga_primer_hits.bed

samtools faidx haemonchus_contortus.PRJEB506.WBPS16.genomic.fa
cat haemonchus_contortus.PRJEB506.WBPS16.genomic.fa.fai | awk '{print $1,$2}' OFS="\t" > haemonchus_contortus.genome

bedtools makewindows -g haemonchus_contortus.genome -w 50000 > haemonchus_contortus.50k_window.bed

bedtools coverage -a haemonchus_contortus.50k_window.bed -b swga_primer_hits.bed > haemonchus_contortus.swga_primer_hits.50k.coverage


>swga_primer_hits.50k.coverage
for i in *.coords; do
     cat ${i} | sort -k1,1 -k2,2n | awk -F "\t" '{print $1,$2,$2+1,$3}' OFS="\t" > swga_primer_hits.bed;
     bedtools coverage -a haemonchus_contortus.50k_window.bed -b swga_primer_hits.bed | awk -v SEQUENCE=${i%.coords} '{print $0,SEQUENCE}' OFS="\t" >> swga_primer_hits.50k.coverage;

done


```



## swga
- alternative approach, based on https://doi.org/10.1093/bioinformatics/btx118
- https://github.com/eclarke/swga

```bash

 module load swga/0.4.4--py27heb12742_2


swga init -f trichuris_trichiura.fa -b sth_minus_ttrichuria.fa

swga count --min_size 8 --max_size 12 --force

swga filter

swga find_sets

swga export sets --order_by score --limit 50 > best_sets.txt

```
