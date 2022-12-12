python /u/project/gandalm/shared/apps/sashimi.py/main.py \
  -e chr7:21516697-21521741 \
  -r /u/project/gandalm/shared/refGenomes/hg19/Gencode/v33/gencode.v33lift37.annotation.sorted.gtf.gz \
  --density sashimipy_density.tsv \
  -o sashimipy.pdf \
  --dpi 300 \
  --width 6 \
  --height 1 \
  --show-junction-num \
  --sites 21521233 \
  --included-junctions chr7:21516927-21521542,chr7:21516927-21521120,chr7:21521302-21521542

# somehow for splice junction coordinates, start is the end of exon+2, end is start of next exon
# figured this out by looking into example code on github and PTBP3 gene exon coordinates GENCODE hg38 v40
# also negative numbers on density plot seem to just be for visualization, putting some archs on the other side (?)
