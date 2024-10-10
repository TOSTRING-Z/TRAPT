
# # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # #     Lisa benchmark dataset  # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # lisa 本地 https://github.com/liulab-dfci/lisa2/tree/master
# 文库下载
# http://cistrome.org/~alynch/data/lisa_data/hg38_1000_2.0.h5
# mv -> [conda_path]/envs/lisa_env/lib/python3.x/site-packages/lisa/data/hg38_1000_2.0.h5
conda activate lisa_env

input_dir=input/Lisa
output_dir=new/result/2.11/lisa
tmp_dir=new/result/2.11
type=up

mkdir -p $output_dir/$type

ls $output_dir/$type/*.txt.lisa.tsv|awk -F "$output_dir/$type/" '{print $2}'|awk -F '.txt.lisa.tsv' '{print $1}' > $tmp_dir/lisa.tr 
ls $input_dir/$type|while read file;do
    tr=$(echo $file|awk -F '.txt' '{print $1}')
    if (grep $tr $tmp_dir/lisa.tr);then
        echo $tr exist!
    else
        input=$input_dir/$type/$file
        output=$output_dir/$type/$file
        ls -lh $input
        lisa oneshot hg38 $input -o $output --save_metadata
    fi
done

# Lisa未识别的TRs
header="Rank	sample_id	factor	cell_line	cell_type	tissue	DNase_p_value	ChIP-seq_p_value	H3K27ac_p_value	summary_p_value"
echo $header > new/result/2.11/lisa/down/ETS1@284_down.txt.lisa.tsv
echo $header > new/result/2.11/lisa/down/REST@271_down.txt.lisa.tsv
echo $header > new/result/2.11/lisa/down/REST@272_down.txt.lisa.tsv
echo $header > new/result/2.11/lisa/down/SOX2@17_down.txt.lisa.tsv
echo $header > new/result/2.11/lisa/up/SOX2@17_up.txt.lisa.tsv

# # # # bart 本地 https://zanglab.github.io/bart/index.htm#install
conda activate bart_env

input_dir=input/Lisa
output_dir=new/result/2.11/bart
tmp_dir=new/result/2.11
type=up
bart=bart2

mkdir -p $output_dir/$type

ls $output_dir/$type/*_bart_results.txt|awk -F "$output_dir/$type/" '{print $2}'|awk -F '_bart_results.txt' '{print $1}' > $tmp_dir/bart.tr
ls $input_dir/$type|while read file;do
    tr=$(echo $file|awk -F '.txt' '{print $1}')
    if (grep $tr $tmp_dir/bart.tr);then
        echo $tr exist!
    else
        input=$input_dir/$type/$file
        output=$output_dir/$type
        ls -lh $input
        $bart geneset -i $input -s hg38 --outdir $output
    fi
done