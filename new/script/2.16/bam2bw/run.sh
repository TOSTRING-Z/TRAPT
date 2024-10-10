wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedClip
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
chmod +x liftOver bigWigToBedGraph bedGraphToBigWig


# bam to bigwig
bamCoverage -b analysis/dupbam/"$fasta_name".sort.bam --ignoreDuplicates \
    --skipNonCoveredRegions \
    --normalizeUsing RPKM \
    --binSize 1 -p max -o analysis/bigwig/$fasta_name.bw
# hg19 to hg38
./bigWigToBedGraph R22025310.bw input.bedGraph
./liftOver -minMatch=0.1 -ends=2 input.bedGraph hg19ToHg38.over.chain.gz output.bedGraph unmapped.bed
# 合并重叠区域
sortBed -i output.bedGraph > output.sorted.bedGraph
bedtools merge -i output.sorted.bedGraph -c 4 -o sum > output.merged.bedGraph

./bedGraphToBigWig output.merged.bedGraph hg38.chrom.sizes R22025310_hg38.bw
