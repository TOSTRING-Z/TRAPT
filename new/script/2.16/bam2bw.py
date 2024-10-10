"""
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedClip
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
chmod +x liftOver bigWigToBedGraph bedGraphToBigWig bedClip
"""

import os
import sys
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--workdir", type=str, default="/data/zgr/data/TRAPT/tool/new/script/2.16/bam2bw/")
parser.add_argument("--genome", type=str, default="hg19",help="hg19/hg38/hg19to38/hg38to19")
parser.add_argument("--bam", type=str, default="new/result/2.16/marks/RB_cas.sort.bam")
parser.add_argument("--output_path", type=str, default="new/result/2.16/bw")
parser.add_argument("--save", type=str, default="false")
parser.add_argument("--force", type=str, default="false")
args = parser.parse_args()

os.environ['PATH'] = os.pathsep.join([args.workdir, os.environ.get('PATH', '')])

name = os.path.basename(args.bam).split(".")[0]
bai = f"{args.bam}.bai"
bw_hg19 = os.path.join(args.output_path,f"{name}.hg19.bw")
bedgraph_hg19 = os.path.join(args.output_path,f"{name}.hg19.bedGraph")
bedgraph_hg38 = os.path.join(args.output_path,f"{name}.hg38.bedGraph")
bedgraph_hg38_sorted = os.path.join(args.output_path,f"{name}.hg38.sorted.bedGraph")
bedgraph_hg38_merged = os.path.join(args.output_path,f"{name}.hg38.merged.bedGraph")
bedgraph_hg19_sorted = os.path.join(args.output_path,f"{name}.hg19.sorted.bedGraph")
bedgraph_hg19_merged = os.path.join(args.output_path,f"{name}.hg19.merged.bedGraph")
bedgraph_hg19_filtered = os.path.join(args.output_path,f"{name}.hg19.filtered.bedGraph")
bedgraph_hg19_clipped = os.path.join(args.output_path,f"{name}.hg19.clipped.bedGraph")
unmapped = os.path.join(args.output_path,f"{name}.unmapped.bed")
bw_hg38 = os.path.join(args.output_path,f"{name}.hg38.bw")

if args.genome == "hg19to38":
    if os.path.exists(bw_hg38) and args.force == "false":
        print(f"file {bw_hg38} final!")
        sys.exit(0)
    os.system(f"""
if [ "{args.force}" = "true" ] || [ ! -e "{bai}" ]; then samtools index {args.bam} {bai}; fi
if [ "{args.force}" = "true" ] || [ ! -e "{bw_hg19}" ]; then bamCoverage -b {args.bam} --ignoreDuplicates --skipNonCoveredRegions --normalizeUsing RPKM --binSize 1 -p max -o {bw_hg19}; fi
if [ "{args.force}" = "true" ] || [ ! -e "{bedgraph_hg19}" ]; then bigWigToBedGraph {bw_hg19} {bedgraph_hg19}; fi
if [ "{args.force}" = "true" ] || [ ! -e "{unmapped}" ]; then liftOver -minMatch=0.1 -ends=2 {bedgraph_hg19} {args.workdir}/hg19ToHg38.over.chain.gz {bedgraph_hg38} {unmapped}; fi
if [ "{args.force}" = "true" ] || [ ! -e "{bedgraph_hg38_sorted}" ]; then sortBed -i {bedgraph_hg38} > {bedgraph_hg38_sorted}; fi
if [ "{args.force}" = "true" ] || [ ! -e "{bedgraph_hg38_merged}" ]; then bedtools merge -d -1 -i {bedgraph_hg38_sorted} -c 4 -o sum > {bedgraph_hg38_merged}; fi
if [ "{args.force}" = "true" ] || [ ! -e "{bw_hg38}" ]; then bedGraphToBigWig {bedgraph_hg38_merged} {args.workdir}/hg38.chrom.sizes {bw_hg38}; fi
if [ "{args.save}" = "false" ]; then rm -f {bw_hg19} {bedgraph_hg19} {bedgraph_hg38} {unmapped} {bedgraph_hg38_sorted} {bedgraph_hg38_merged}; fi
""")
if args.genome == "hg38to19":
    if os.path.exists(bw_hg19) and args.force == "false":
        print(f"file {bw_hg19} final!")
        sys.exit(0)
    os.system(f"""
if [ "{args.force}" = "true" ] || [ ! -e "{bai}" ]; then samtools index {args.bam} {bai}; fi
if [ "{args.force}" = "true" ] || [ ! -e "{bw_hg38}" ]; then bamCoverage -b {args.bam} --ignoreDuplicates --skipNonCoveredRegions --normalizeUsing RPKM --binSize 1 -p max -o {bw_hg38}; fi
if [ "{args.force}" = "true" ] || [ ! -e "{bedgraph_hg38}" ]; then bigWigToBedGraph {bw_hg38} {bedgraph_hg38}; fi
if [ "{args.force}" = "true" ] || [ ! -e "{unmapped}" ]; then liftOver -minMatch=0.1 -ends=2 {bedgraph_hg38} {args.workdir}/hg38ToHg19.over.chain.gz {bedgraph_hg19} {unmapped}; fi
if [ "{args.force}" = "true" ] || [ ! -e "{bedgraph_hg19_sorted}" ]; then sortBed -i {bedgraph_hg19} > {bedgraph_hg19_sorted}; fi
if [ "{args.force}" = "true" ] || [ ! -e "{bedgraph_hg19_merged}" ]; then bedtools merge -d -1 -i {bedgraph_hg19_sorted} -c 4 -o sum > {bedgraph_hg19_merged}; fi
if [ "{args.force}" = "true" ] || [ ! -e "{bedgraph_hg19_filtered}" ]; then chromosome_filter {bedgraph_hg19_merged} > {bedgraph_hg19_filtered}; fi
if [ "{args.force}" = "true" ] || [ ! -e "{bedgraph_hg19_clipped}" ]; then bedClip {bedgraph_hg19_filtered} {args.workdir}/hg19.chrom.sizes {bedgraph_hg19_clipped}; fi
if [ "{args.force}" = "true" ] || [ ! -e "{bw_hg19}" ]; then bedGraphToBigWig {bedgraph_hg19_clipped} {args.workdir}/hg19.chrom.sizes {bw_hg19}; fi
if [ "{args.save}" = "false" ]; then rm -f {bw_hg38} {bedgraph_hg38} {bedgraph_hg19} {unmapped} {bedgraph_hg19_sorted} {bedgraph_hg19_merged} {bedgraph_hg19_filtered} {bedgraph_hg19_clipped}; fi
""")
if args.genome == "hg38":
    if os.path.exists(bw_hg38) and args.force == "false":
        print(f"file {bw_hg38} final!")
        sys.exit(0)
    os.system(f"bamCoverage -b {args.bam} --ignoreDuplicates --skipNonCoveredRegions --normalizeUsing RPKM --binSize 1 -p max -o {bw_hg38}")
if args.genome == "hg19":
    if os.path.exists(bw_hg19) and args.force == "false":
        print(f"file {bw_hg19} final!")
        sys.exit(0)
    os.system(f"bamCoverage -b {args.bam} --ignoreDuplicates --skipNonCoveredRegions --normalizeUsing RPKM --binSize 1 -p max -o {bw_hg19}")
    