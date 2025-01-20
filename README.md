## TRAPT

[TRAPT](https://bio.liclab.net/TRAPT) is a multi-stage fused deep learning framework for transcription regulators prediction via integraing large-scale epigenomic data.

## Usage

### Online
We have developed a corresponding web service (https://bio.liclab.net/TRAPT). The website is designed to accepts gene sets input by users for analysis, allowing easy retrieval of analytical results. We've also thoughtfully included an email notification feature. On the results page, the website displays all TR activity scores, as well as the ranking and all individual scores of TRs. Concurrently, the website provides annotation details and relevant quality control information for each transcriptional regulator. Compared to offline tools, online analysis tools offer additional features on the browsing and result analysis pages. The online tools facilitate visualization of the predicted 3D protein structure for each TR, leveraging AlphaFold's predictions. Additionally, the online tools incorporate a genome browser to facilitate user interaction with the genomic tracks associated with each TR.

### Offline

First, download library: 

[TRAPT library](https://bio.liclab.net/TRAPT/download) or 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8080171.svg)](https://doi.org/10.5281/zenodo.8080171)




Second, install TRAPT:

```shell
conda create --name TRAPT python=3.7
conda activate TRAPT
pip install TRAPT
```

Run TRAPT using a [case](https://bio.liclab.net/TRAPT/static/download/ESR1@DataSet_01_111_down500.txt):

Use shell commands:
```shell
# help
trapt --help
# run
trapt --library library \
      --input ESR1@DataSet_01_111_down500.txt \
      --output output/test/ESR1@DataSet_01_111_down500
```

Using the python interface:
```python
import os
from TRAPT.Tools import Args, RP_Matrix
from TRAPT.Run import runTRAPT

# library path
library = 'library'
# input file path
input = 'ESR1@DataSet_01_111_down500.txt'
# output file path
output = 'output/test/ESR1@DataSet_01_111_down500'

args = Args(input, output, library)
runTRAPT(args)
```

### De novo library

```bash
# Reconstruct TR-H3K27ac adjacency matrix
python3 src/TRAPT/DLVGAE.py H3K27ac library
# Reconstruct TR-ATAC adjacency matrix
python3 src/TRAPT/DLVGAE.py ATAC library
# Prediction D-RP(H3K27ac) matrix
python3 src/TRAPT/CalcTRSampleRPMatrix.py H3K27ac library
# Prediction D-RP(ATAC) matrix
python3 src/TRAPT/CalcTRSampleRPMatrix.py ATAC library
```

### Parameter Documentation
Documentation can be found on [readthedocs](https://trapt.readthedocs.io/en/latest/) or [pdf](https://github.com/TOSTRING-Z/TRAPT/blob/main/docs/_build/latex/trapt.pdf).

## Case 1: Inferring Key Regulators in Breast Cancer Using Bulk RNA-seq
### Data Preparation
```shell
mkdir -p new/result/case1
wget -c https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz -O new/result/case1/TCGA.BRCA.sampleMap_HiSeqV2.gz
```
### Differential Expression Analysis
```shell
Rscript new/script/case1/data_processing.r
```

<img src="./new/result/case1/volcano_plot_limma_voom.svg" alt="" style="width: 50%" />

### TRAPT Infers Key Regulators
```shell
trapt --library library \
      --input new/result/case1/genes.txt \
      --output new/result/case1/trapt_output
```

## Case 2: Inferring key regulators in human hematopoietic stem cells using scRNA-seq data.
### Data Preparation
```shell
mkdir -p new/result/case2
wget -c https://bio.liclab.net/TRAPT/static/download/data.csv.gz -O new/result/case2/data.csv.gz
```
### Data Processing (Dimensionality Reduction, Clustering, Identification of Marker Genes)
```shell
Rscript new/script/case2/data_processing.r
python3 new/script/case2/integrated_degs.py
Rscript new/script/case2/multi_volcano_plot.r
```

<img src="./new/result/case2/cluster-umap.svg" alt="" style="width: 50%;" />

<img src="./new/result/case2/cluster-pca.svg" alt="" style="width: 50%" />

<img src="./new/result/case2/mult_volcano_plot.svg" alt="" style="width: 50%" />


### TRAPT Infers Key Regulators of Cell Fate Determination
```shell
trapt --library library \
      --input 'new/result/case2/output/markers_top200-[pDCs-LMPs].txt' \
      --output new/result/case2/trapt_output
```


## Research details

### The differences in model variant predictions TR
```shell
python3 new/script/2.8/RankTRAPT_var.py --output_dir output/KnockTFv1 --rank_path new/result/2.8/files --name TRAPT-H3K27ac --model H3K27ac
python3 new/script/2.8/RankTRAPT_var.py --output_dir output/KnockTFv1 --rank_path new/result/2.8/files --name TRAPT-ATAC --model ATAC

python3 new/script/2.8/var_predict.py --rank_x new/result/2.8/files/rank_TRAPT-H3K27ac.csv --rank_y new/result/2.8/files/rank_TRAPT-ATAC.csv --output_path new/result/2.8/figure/rank_scatterplot_h3k27ac-atac.svg --title "TRAPT-H3K27ac vs. TRAPT-ATAC"
```

<img src="./new/result/2.8/figure/rank_scatterplot_h3k27ac-atac.svg" alt="" style="width: 50%" />


```shell
python3 new/script/2.8/RankTRAPT_var.py --output_dir output/KnockTFv1 --rank_path new/result/2.8/files --name TRAPT --model TRAPT
python new/script/2.8/TR_RP_run.py --input_path input/KnockTFv1/down --output_path new/result/2.8/output-TR_RP
python new/script/2.8/TR_RP_run.py --input_path input/KnockTFv1/up --output_path new/result/2.8/output-TR_RP
python3 new/script/2.8/RankTRAPT_var.py --output_dir new/result/2.8/output-TR_RP --rank_path new/result/2.8/files --name TR-model --model TR
# peak in promoter
python3 new/script/2.8/peak_in_promoter_matrix.py --library library --output library

python new/script/2.8/TR_Peak_run.py --input_path input/KnockTFv1/down --output_path new/result/2.8/output-TR_peak
python new/script/2.8/TR_Peak_run.py --input_path input/KnockTFv1/up --output_path new/result/2.8/output-TR_peak
python3 new/script/2.8/RankTRAPT_var.py --output_dir new/result/2.8/output-TR_peak --rank_path new/result/2.8/files --name Peak-promoter-model --model Peak

python3 new/script/2.8/var_predict.py --rank_x new/result/2.8/files/rank_TRAPT.csv --rank_y new/result/2.8/files/rank_Peak-promoter-model.csv --output_path new/result/2.8/figure/rank_scatterplot_trapt-peak.svg --title "TRAPT vs. Peak-promoter-model"
```

<img src="./new/result/2.8/figure/rank_scatterplot_trapt-peak.svg" alt="" style="width: 50%" />

```shell
python3 new/script/2.8/Rank_var.py --output_path new/result/2.8/figure --name Model-variant --type ALL --rank_path new/result/2.8/files --source KnockTF --columns TRAPT,TRAPT-H3K27ac,TR-model,Peak-promoter-model --col_names 'TRAPT,TRAPT(ATAC-),TRAPT(H3K27ac-),Peak-promoter-model'
```

<img src="./new/result/2.8/figure/rank_Model-variant@mmr_bar.svg" alt="" style="width: 50%" />


### Algorithm comparison using cistrome background data
```shell
# KnockTF benchmark dataset
python3 new/script/2.11/RankTRAPT_Source.py --output_dir output/KnockTFv1 --rank_path new/result/2.11/files --source cistrome --name TRAPT-source_of_cistrome --model TRAPT
## Lisa
python3 new/script/2.11/RankLisa.py --match_dir output/KnockTFv1 --input_path other/lisa --type down,up --output_path new/result/2.11/files
## BART
python3 new/script/2.11/RankBART.py --match_dir output/KnockTFv1 --input_path other/bart --type down,up --output_path new/result/2.11/files

python3 new/script/3.11/Rank_all.py --output_path new/result/2.11/figure --name TRAPT-cistrome_TF --type TF --rank_path new/result/2.11/files --source KnockTF --columns TRAPT-source_of_cistrome,Lisa,BART --col_names 'TRAPT-source_of_cistrome,Lisa,BART'

python3 new/script/3.11/Rank_all.py --output_path new/result/2.11/figure --name TRAPT-cistrome_CR_TcoF --type CR_TcoF --rank_path new/result/2.11/files --source KnockTF --columns TRAPT-source_of_cistrome,Lisa,BART --col_names 'TRAPT-source_of_cistrome,Lisa,BART'

python3 new/script/2.11/Rank_pile.py --output_path new/result/2.11/figure --name TRAPT-cistrome --type ALL --rank_path new/result/2.11/files --source KnockTF --columns TRAPT-source_of_cistrome,Lisa,BART --col_names 'TRAPT-source_of_cistrome,Lisa,BART'
```

<img src="./new/result/2.11/figure/rank_TRAPT-cistrome@boxplot.svg" alt="" style="width: 50%" />

<img src="./new/result/2.11/figure/rank_TRAPT-cistrome@mmr_bar.svg" alt="" style="width: 50%" />

<img src="./new/result/2.11/figure/rank_TRAPT-cistrome@mmr_pile_bar.svg" alt="" style="width: 50%" />

```shell
# Lisa benchmark dataset

## TRAPT
python3 new/script/2.11/run_mult.py --input_dir input/Lisa/down --output_dir new/result/2.11/output-Lisa
python3 new/script/2.11/run_mult.py --input_dir input/Lisa/up --output_dir new/result/2.11/output-Lisa
python3 new/script/2.11/RankTRAPT_Source.py --output_dir new/result/2.11/output-Lisa --rank_path new/result/2.11/files_lisa --source cistrome --name TRAPT-cistrome_lisa --model TRAPT

## Lisa new/script/2.11/lisa_bart.sh
python3 new/script/2.11/RankLisa.py --match_dir new/result/2.11/output-Lisa --input_path new/result/2.11/lisa --type down,up --output_path new/result/2.11/files_lisa

## BART new/script/2.11/lisa_bart.sh
python3 new/script/2.11/RankBART.py --match_dir new/result/2.11/output-Lisa --input_path new/result/2.11/bart --type down,up --output_path new/result/2.11/files_lisa

## ChEA3
python3 new/script/2.11/chea3_run.py --input_path input/Lisa/down --output_path new/result/2.11/chea3/down
python3 new/script/2.11/chea3_run.py --input_path input/Lisa/up --output_path new/result/2.11/chea3/up
python3 new/script/2.11/RankChEA3.py --match_dir new/result/2.11/output-Lisa --input_path new/result/2.11/chea3 --type down,up --output_path new/result/2.11/files_lisa

## i-cisTarget
python3 new/script/2.11/icistarget_run.py --input_path input/Lisa/down --output_path new/result/2.11/icistarget/down --download_dir /root/下载 --executable_path /install/chromedriver_linux64/chromedriver
python3 new/script/2.11/icistarget_run.py --input_path input/Lisa/up --output_path new/result/2.11/icistarget/up --download_dir /root/下载 --executable_path /install/chromedriver_linux64/chromedriver
python3 new/script/2.11/RankicisTarget.py --match_dir new/result/2.11/output-Lisa --input_path new/result/2.11/icistarget --type down,up --output_path new/result/2.11/files_lisa

python3 new/script/2.11/Rank.py --output_path new/result/2.11/figure --name Lisa-benchmark-dataset --type ALL --rank_path new/result/2.11/files_lisa --source Lisa --columns TRAPT-cistrome_lisa,Lisa,BART,i-cisTarget,ChEA3
```

<img src="./new/result/2.11/figure/rank_Lisa-benchmark-dataset@boxplot.svg" alt="" style="width: 50%" />

<img src="./new/result/2.11/figure/rank_Lisa-benchmark-dataset@mmr_bar.svg" alt="" style="width: 50%" />


### TR decay rate
```shell
python3 new/script/2.5/CalcTRRPMatrix.py --library library --output new/result/2.5
Rscript new/script/2.5/gsea_run.r
python3 new/script/2.5/decay_rate.py
```

<img src="./new/result/2.5/long.svg" alt="" style="width: 50%" />

<img src="./new/result/2.5/short.svg" alt="" style="width: 50%" />


### Overlap between selected epigenetic profiles and TR peaks
```shell
sed '1d' input/AD/genename_out.txt | awk '{print $10}' | head -n 200 > new/result/2.16/input/AD_200.txt
python3 src/TRAPT/Run.py --input new/result/2.16/input/AD_200.txt --output new/result/2.16/output/AD_200
# 172.27.0.239 /data/zgr/tmp
python3 get_marks.py --ssh_user zhangguorui --ssh_ip 10.100.3.6 --marks_path /wyzdata8/SEdb1/all_bam/sample --input H3K27ac_info_samples.csv --output_path marks
python3 get_marks.py --ssh_user zhangguorui --ssh_ip 10.100.3.6 --marks_path /wyzdata5/SEdb2/bam --input H3K27ac_info_samples.csv --output_path marks

python3 new/script/2.16/bam2bw.py
ls new/result/2.16/marks/hg19/*.bam|while read bam;do
    python3 new/script/2.16/bam2bw.py --bam $bam
done
ls new/result/2.16/marks/hg38/*.bam|while read bam;do
    python3 new/script/2.16/bam2bw.py --genome hg38to19 --bam $bam
done
# IGV visualization
```

### Dataset analysis code, chart code
```shell
# Fig. 2 | Evaluation of TRAPT and comparative methods on TR knockdown/knockout and TF binding datasets.
### TRAPT
python3 new/script/2.11/RankTRAPT_Source.py --output_dir output/KnockTFv1 --name TRAPT --type down500,up500 --rank_path 'new/result/3.11/Fig. 2/files'
### Lisa
python3 new/script/2.11/RankLisa.py --match_dir output/KnockTFv1 --input_path other/lisa --type down,up --output_path 'new/result/3.11/Fig. 2/files'
### BART
python3 new/script/2.11/RankBART.py --match_dir output/KnockTFv1 --input_path other/bart --type down,up --output_path 'new/result/3.11/Fig. 2/files'
### ChEA3
python3 new/script/2.11/RankChEA3.py --match_dir output/KnockTFv1 --input_path other/chea3 --type down,up --output_path 'new/result/3.11/Fig. 2/files'
### i-cisTarget
python3 new/script/2.11/RankicisTarget.py --match_dir output/KnockTFv1 --input_path other/icistarget --type down,up --output_path 'new/result/3.11/Fig. 2/files'
### Plotting
python3 new/script/3.11/Rank_all.py --output_path 'new/result/3.11/Fig. 2/figure' --name rank_KnockTF-benchmark-dataset --type ALL --rank_path 'new/result/3.11/Fig. 2/files' --source KnockTF --columns TRAPT,Lisa,BART,i-cisTarget,ChEA3
```

<img src="./new/result/3.11/Fig.%202/figure/rank_rank_KnockTF-benchmark-dataset@mmr_bar.svg" alt="" style="width: 50%" />

<img src="./new/result/3.11/Fig.%202/figure/rank_rank_KnockTF-benchmark-dataset@boxplot.svg" alt="" style="width: 50%" />

<img src="./new/result/3.11/Fig.%202/figure/rank_rank_KnockTF-benchmark-dataset@bar.svg" alt="" style="width: 50%" />

<img src="./new/result/3.11/Fig.%202/figure/rank_rank_KnockTF-benchmark-dataset@line.svg" alt="" style="width: 50%" />

<img src="./new/result/3.11/Fig.%202/figure/rank_rank_KnockTF-benchmark-dataset@groupbar.svg" alt="" style="width: 50%" />

```shell
# Fig. 3 | Using the differential gene sets from TR knockdown/knockout experiments by KnockTF, we evaluated the performance of TRAPT.
python3 new/script/3.11/rank_scatter_plot.py
```

<img src="./new/result/3.11/Fig.%203/figure/rank_scatter.svg" alt="" style="width: 50%" />

```shell
# Fig. 4 | Illustration of the TRAPT framework using downregulated genes from ESR1 gene knockout experiments in gastric cancer and KMCF7 breast cancer.
python3 new/script/3.11/case_esr1.py
```

<img src="./new/result/3.11/Fig.%204/figure/ESR1-KD-MCF7.svg" alt="" style="width: 50%" />

<img src="./new/result/3.11/Fig.%204/figure/ppi.svg" alt="" style="width: 50%" />

<img src="./new/result/3.11/Fig.%204/GTEx.svg" alt="" style="width: 50%" />

<img src="./new/result/3.11/Fig.%204/TCGA.svg" alt="" style="width: 50%" />

```shell
# Fig. 5 | Prediction of functional transcriptional regulators for Alzheimer's disease using post-GWAS analysis.
python3 new/script/3.11/enrichment.py
```

<img src="./new/result/3.11/Fig.%205/figure/enrichment.svg" alt="" style="width: 50%" />


```shell
# Fig. 6 | TRAPT identifies transcriptional regulators associated with cell fate and tissue identity. 
```shell
python3 new/script/case2/gtex_diff_analysis.py
python3 new/script/3.11/gtex_heatmap.py
```

<img src="./new/result/3.11/Fig.%206/figure/GTEx-heatmap.svg" alt="" style="width: 50%" />
