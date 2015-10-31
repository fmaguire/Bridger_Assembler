#!/bin/bash -ve


if [ -e reads.right.fq.gz ] && [ ! -e reads.right.fq ]; then
    gunzip -c reads.right.fq.gz > reads.right.fq
fi

if [ -e reads.left.fq.gz ] && [ ! -e reads.left.fq ]; then
    gunzip -c reads.left.fq.gz > reads.left.fq
fi


if [! -e bridger_out_dir/Bridger.fasta]; then
    ./run_Me_as_DS.sh
fi


# use RSEM to estimate read abundance

../perllib/run_RSEM_align_n_estimate.pl  --transcripts bridger_out_dir/Bridger.fasta --seqType fq --left reads.left.fq --right reads.right.fq
