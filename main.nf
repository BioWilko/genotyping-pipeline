data_dir = "" //directory containing artic pipeline output
out_dir = ""
metadata = "" //path to TSV containing all relevant metadata formatted for CLIMB-COVID
definition_dir = "" //directory containing aln2type yaml files, e.g. ./definitions in this repo

samp_ids = params.id_list
id_str = samp_ids.replaceAll(",", " ")

process custom_metadata {
  output:
    file("custom_metadata.tsv") into metadata_ch

  """
  id_array=(${id_str})
  head -1 ${metadata} > custom_metadata.tsv
  for i in "\${id_array[@]}"
  do
    grep "\$i" ${metadata} >> custom_metadata.tsv || true
  done
  """
}

process pull_cog_ids {
  input:
    file metadata from metadata_ch
  output:
    file("cog_ids.csv") into id_ch
    file("sender_ids.tsv") into sender_ids_ch

  """
  #!/usr/bin/env python
  import pandas as pd
  import numpy as np

  metadata_df = pd.read_csv("${metadata}", sep="\t")
  metadata_df = metadata_df.drop_duplicates(subset="central_sample_id")
  metadata_df = metadata_df.loc[metadata_df["exclude"].isnull()]
  metadata_df["coverage"] = ""
  id_df = metadata_df.copy()

  for index, row in metadata_df.iterrows():
    try:
      coverage_df = pd.read_csv("${data_dir}" + row["central_sample_id"] + ".cov.txt", index_col="sample", sep="\t")
      metadata_df.at[index, "coverage"] = coverage_df.at[(row["central_sample_id"]),"perc"]
    except:
      id_df.drop(labels=index, inplace=True)
      metadata_df.at[index, "coverage"] = "N/A"


  id_df["central_sample_id"].to_csv("cog_ids.csv", index=False, header=False)
  metadata_df[["central_sample_id", "sender_sample_id", "coverage"]].to_csv("sender_ids.tsv", index=False, sep="\t")
  """
}

id_ch
  .splitText()
  .map{it -> it.trim()}
  .into {aln_ids_ch;consensus_id_ch}

process cat_msa {
  input:
    val samp_ids from aln_ids_ch.toList()
  output:
    file("combined_msa.fasta") into msa_out_ch
  script:
    id_str = samp_ids.join(',')

  """
  cat ${data_dir}{${id_str}}.nanopolish-indel.muscle.out.fasta > temp_msa.fasta
  awk -F'.' '/^>/ {print \$1; next}{print}' < temp_msa.fasta > combined_msa.fasta
  """    
}

process cat_consensus {
  input:
    val samp_ids from consensus_id_ch.toList()
  output:
    file("combined_consensus.fasta") into consensus_out_ch
  script:
    id_str = samp_ids.join(',')

  """
  cat ${data_dir}{${id_str}}.nanopolish-indel.consensus.fasta > temp_consensus.fasta
  awk -F'.' '/^>/ {print \$1; next}{print}' < temp_consensus.fasta > combined_consensus.fasta
  """
}

process aln2type {
  input:
    file msa from msa_out_ch
  output:
    file("type_calls.csv") into type_ch
  
  """
  aln2type ./ ./ type_calls.csv MN908947 ${msa} --no_call_deletion --output_unclassified /data/homes/samw/projects/turnkey_genotype_validation/definitions/*
  """
}

process pangolin {
  conda '/data/homes/samw/miniconda3/envs/pangolin/'
  input:
    file combined from consensus_out_ch
  output:
    file("lineage_report.csv") into lineage_ch

  """
  pangolin ${combined}
  """
}

process report {
  publishDir out_dir, overwrite: true
  input:
    file types from type_ch
    file lineage from lineage_ch
    file sender_ids from sender_ids_ch
  output:
    file("turnkey_genotype_calls.tsv")

  """
  #!/usr/bin/env python
  import pandas as pd
  from Bio import AlignIO
  import glob
  import os

  df = pd.read_csv("${types}")

  group_df = df.groupby(["sample_id"])

  ids = df["sample_id"].unique()

  # mutations = ["E484K", "E484Q", "N501Y", "K417N", "K417T", "L452R", "P681R", "P681H"]

  dictionary = {}

  mutation_locs = {}

  definitions = glob.glob("${definition_dir}/*.yml")

  mutations = [x.split("/")[-1].replace(".yml", "") for x in definitions]

  for definition_file in definitions:
    steam = os.popen('grep "one-based-reference-position:" %s' %(definition_file))
    location = steam.read().split(": ")[1]
    mutation = definition_file.split("/")[-1].replace(".yml", "")
    mutation_locs[mutation] = int(location)

  for samp_id in ids:
    samp_df = group_df.get_group(samp_id)
    alignment = AlignIO.read(open("${data_dir}/%s.nanopolish-indel.muscle.out.fasta" %(samp_id)), "fasta")
    if samp_df["phe-label"].any() != "unclassified":
      dictionary[samp_id] = { i : "" for i in mutations }
      for mutation in mutations:
          if alignment[0].seq[mutation_locs[mutation]] != "N":
            dictionary[samp_id][mutation] = "N"
          else:
            dictionary[samp_id][mutation] = "X"
      for index, row in samp_df.iterrows():
        if row["status"] == "confirmed":
          dictionary[samp_id][row["phe-label"]] = "Y"

  unclassified = df.loc[df["phe-label"] == "unclassified"]

  for index, row in unclassified.iterrows():
    dictionary[row["sample_id"]] = dict.fromkeys(mutations, "N")

  type_calls = pd.DataFrame.from_dict(dictionary, orient='index')
  lineage_df = pd.read_csv("${lineage}", index_col="taxon")
  lineage_df = lineage_df["lineage"]
  sender_ids = pd.read_csv("${sender_ids}", index_col="central_sample_id", sep="\t")

  lineage_df = sender_ids.join(lineage_df, how="outer")

  report_df = lineage_df.join(type_calls, how="outer")

  report_df.to_csv("turnkey_genotype_calls.tsv", sep="\t")
  """
}


