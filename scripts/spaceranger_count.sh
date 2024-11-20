PREFIX=/home/ubuntu/mnt/bhoomi-bcl-dl/
SAMP=Slide_1
CYTAIMAGE=${PREFIX}/images/CAVG10138R1_2024-08-29_12-06-09_lynch_bhb-visiumhd-8-29_H1-J4XFQXR_A1_slide1.tif
IMAGE=${PREFIX}/images/Slide_1.tif
SLIDE=H1-J4XFQXR
AREA=A1
CBS=64

# SAMP=Slide_2
# CYTAIMAGE=${PREFIX}/images/CAVG10138R1_2024-08-29_12-06-09_lynch_bhb-visiumhd-8-29_H1-J4XFQXR_D1_slide2.tif
# IMAGE=${PREFIX}/images/Slide_2.tif
# SLIDE=H1-J4XFQXR
# AREA=D1

~/pkgs/spaceranger-3.1.1/spaceranger count \
    --fastqs=${PREFIX}/fastq \
    --id=bhoomi_bhb_repeat_${SAMP}_HE_bs${CBS} \
    --sample=${SAMP} \
    --output-dir=${PREFIX}/spaceranger/${SAMP}_HE_bs${CBS} \
    --cytaimage=${CYTAIMAGE} \
    --image=${IMAGE} \
    --slide=${SLIDE}\
    --area=${AREA} \
    --probe-set=${PREFIX}/ref/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
    --transcriptome=${PREFIX}/ref/refdata-gex-GRCh38-2020-A \
    --create-bam=true \
    --localcores=15 \
    --localmem=100 \
    --custom-bin-size $CBS
