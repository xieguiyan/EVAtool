rm -r /home/liucj/tmp/DRR006758


/home/liucj/tools/anaconda3/bin/python \
  /home/liucj/github/pipelines/pipelines/wdl-smRNA/miRNA_scripts/mir_pipeline_zhangq_20180902.py \
  -i /home/liucj/tmp/DRR006758.sra \
  -n DRR006758 \
  -s /home/liucj/tmp/DRR006758.sh \
  -c /home/liucj/github/pipelines/pipelines/wdl-smRNA/miRNA_scripts/config.miRNA.s10.txt \
  -d /home/liucj/tmp/DRR006758

bash /home/liucj/tmp/DRR006758.sh