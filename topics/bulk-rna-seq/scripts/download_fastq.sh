#!/usr/bin/env bash
# Download all 18 RNA-seq runs of GSE182190 and rename to sample names
set -uo pipefail
OUT=$1   # e.g. data/fastq_files (run from the topic folder, log to data/logs/download.log)
mkdir -p "$OUT" "$OUT/sra"
cat > "$OUT/samples.tsv" <<'S'
SRR15477578	FOXA1_U1
SRR15477580	FOXA1_U2
SRR15477582	FOXA1_U3
SRR15477584	FOXA1_I1
SRR15477586	FOXA1_I2
SRR15477588	FOXA1_I3
SRR15477590	HNF4A_U1
SRR15477592	HNF4A_U2
SRR15477594	HNF4A_U3
SRR15477596	HNF4A_I1
SRR15477598	HNF4A_I2
SRR15477600	HNF4A_I3
SRR15477602	FOXA1HNF4A_U1
SRR15477604	FOXA1HNF4A_U2
SRR15477606	FOXA1HNF4A_U3
SRR15477608	FOXA1HNF4A_I1
SRR15477611	FOXA1HNF4A_I2
SRR15477613	FOXA1HNF4A_I3
S
dl() {
  srr=$1; name=$2; out=$3
  [ -s "$out/${name}_2.fastq.gz" ] && { echo "skip $name"; return; }
  for try in 1 2 3; do
    prefetch "$srr" -O "$out/sra" --max-size 20g \
      && fasterq-dump "$out/sra/$srr" --split-files -e 4 -O "$out" -t "$out/tmp_$srr" \
      && break
    echo "retry $srr ($try)"
  done
  pigz -p 4 -c "$out/${srr}_1.fastq" > "$out/${name}_1.fastq.gz" && rm "$out/${srr}_1.fastq"
  pigz -p 4 -c "$out/${srr}_2.fastq" > "$out/${name}_2.fastq.gz" && rm "$out/${srr}_2.fastq"
  rm -rf "$out/tmp_$srr" "$out/sra/$srr"
  echo "done $name"
}
export -f dl
cut -f1,2 "$OUT/samples.tsv" | xargs -P 3 -n 2 bash -c 'dl "$0" "$1" "'"$OUT"'"'
