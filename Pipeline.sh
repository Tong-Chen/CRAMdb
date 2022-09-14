[TOC]
### database=db and work directory=wd
   
    db=/c/public
    wd=/e/PRJNA222384
    PATH=$PATH:${db}/win
    cd ${wd}


    mkdir -p temp

    cat -A result/metadata.txt | head -n3
   
    sed -i 's/\r/\n/' result/metadata.txt
    cat -A result/metadata.txt | head -n3

    gunzip seq/*.gz
    
### 1. Merge paired reads and rename sequences
   
    for i in `tail -n+2 result/metadata.txt | cut -f 1`;do
      vsearch --fastq_mergepairs seq/${i}_1.fastq. --reverse seq/${i}_2.fastq. \
      --fastqout temp/${i}.merged.fq --relabel ${i}.
    done &

    # If it's not paired reads
    # for i in `tail -n+2 result/metadata.txt | cut -f 1`;do
    # usearch -fastx_relabel seq/${i}.fastq -fastqout temp/${i}.merged.fq -prefix ${i}.
    # done &
  
    cat temp/*.merged.fq > temp/all.fq
   
    less temp/all.fq


### 2.  Quality control (Cut primers and quality filter)

    time vsearch --fastx_filter temp/all.fq \
      --fastq_stripleft 0 --fastq_stripright 0 \
      --fastq_maxee_rate 0.01 \
      --fastaout temp/filtered.fa
      

    vsearch --derep_fulllength temp/filtered.fa \
      --minuniquesize 8 --sizeout --relabel Uni_ \
      --output temp/uniques.fa 

### 3  Denoise: predict biological sequences and filter chimeras
   
    time usearch -unoise3 temp/uniques.fa \
      -zotus temp/zotus.fa
   
    sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa

    mkdir -p result/raw

    cp -f temp/otus.fa result/raw/otus.fa
   
    sed -i 's/\r//g' result/raw/otus.fa


### 4.Generating a feature table

    time vsearch --usearch_global temp/filtered.fa \
      --db result/raw/otus.fa \
      --id 0.97 --threads 4 \
    	--otutabout result/raw/otutab.txt 
    
    sed -i 's/\r//' result/raw/otutab.txt

### 5.Remove plastids and non-bacteria

    # RDP:(rdp_16s_v16_sp) UNITE:(utax_reference_dataset_04.02.2020.fasta.gz)
    
    vsearch --sintax result/raw/otus.fa \
      --db ${db}/usearch/utax_reference_dataset_04.02.2020.fasta.gz \
      --sintax_cutoff 0.6 \
      --tabbedout result/raw/otus.sintax 

    cp result/raw/otu* result/

### 5.1 Summary OTUs table
    usearch -otutab_stats result/otutab.txt \
      -output result/otutab.stat
    cat result/otutab.stat

### 6.Normlize by subsample

    mkdir -p result/alpha
    Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
      --depth 0 --seed 1 \
      --normalize result/otutab_rare.txt \
      --output result/alpha/vegan.txt
    usearch -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat
    cat result/otutab_rare.stat

### 7.Taxonomic summary of species annotations

    cut -f 1,4 result/otus.sintax \
      |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
      > result/taxonomy2.txt

    
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
      result/taxonomy2.txt > temp/otus.tax
    sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
      sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/taxonomy.txt


    mkdir -p result/tax
    sed -i 's/\"//g' result/otus.sintax
    for i in p c o f g s;do
      usearch -sintax_summary result/otus.sintax \
      -otutabin result/otutab_rare.txt -rank ${i} \
      -output result/tax/sum_${i}.txt
    done
    sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
    

### 8.Difference comparison

    mkdir -p result/compare/
    
    compare="Health-Diarrhea"
    Rscript ${db}/script/compare.R \
      --input result/otutab.txt --design result/metadata.txt \
      --group Group --compare ${compare} --threshold 0.1 \
      --method edgeR --pvalue 0.05 --fdr 0.2 \
      --output result/compare/

