# genotyping-pipeline

Requires [aln2type](https://github.com/connor-lab/aln2type), note that this pipeline is provided for transparancy only and has not been tested on other systems, use at your own risk.

Also requires a working install of [pangolin](https://github.com/cov-lineages/pangolin), point towards your conda environment with the flag "--pango".

Unfortunately this will fall over if you only provide one findable COG-ID.

## Usage Example

```
nextflow main.nf --id_list TEST-123456,TEST-123457 --data ~/artic_out_dir/ --out ~/out_dir/ --metadata ~/metadata.tsv --pango ~/miniconda3/envs/pangolin
```