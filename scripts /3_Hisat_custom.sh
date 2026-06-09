#!/usr/bin/env bash
#####################################
# Usage & parsing
#####################################
usage() {

  cat <<'EOF'

HiSat mapping for hints generation 28_8_2025 (f.fontana)
Questo script mappa dati RNA-seq sul genoma e genera junctions/hints per AUGUSTUS.

USO:

  ./hisat2_hints.sh -x <genome.fa> [-s <idx_prefix>] [-t <threads>] [-m <max_intron>] [-l <lista_SE>] [-L <lista_PE>]

Opzioni:
  -x  (obbligatorio)  FASTA del genoma
  -s  (opzionale)     prefisso index hisat2 (default: idx_Genome)
  -t  (opzionale)     threads (default: 40)
  -m  (opzionale)     max intron length (default: 100000)
  -l  (opzionale)     file lista Single-End; una path per riga
  -L  (opzionale)     file lista Paired-End; per riga: <R1><TAB><R2>
  -h                  help

Note:
  * È sufficiente una delle due liste (-l o -L).
  * La lista PE DEVE avere R1 e R2 separati da TAB in ogni riga.

EOF

}

# defaults
idx_Genome_Name="idx_Genome"
genome=""
Threads=40
MaxIntronLength=100000
listsingle=""
listpaired=""

# parsing
while getopts ":x:s:t:m:l:L:h" opt; do

  case "$opt" in
    x) genome="$OPTARG" ;;
    s) idx_Genome_Name="$OPTARG" ;;
    t) Threads="$OPTARG" ;;
    m) MaxIntronLength="$OPTARG" ;;
    l) listsingle="$OPTARG" ;;
    L) listpaired="$OPTARG" ;;
    h) usage; exit 0 ;;
    \?) echo "Errore: opzione non valida -$OPTARG" >&2; usage; exit 2 ;;
    :)  echo "Errore: l'opzione -$OPTARG richiede un argomento." >&2; usage; exit 2 ;;

  esac

done

shift $((OPTIND-1))


# controlli input
[[ -n "$genome" ]] || { echo "Errore: -x <genome.fa> è obbligatorio." >&2; usage; exit 2; }
[[ -s "$genome" ]] || { echo "Errore: file FASTA non trovato o vuoto: $genome" >&2; exit 2; }

if [[ -z "$listsingle" && -z "$listpaired" ]]; then
  echo "Errore: serve almeno una lista (-l o -L)." >&2; usage; exit 2
fi

# validazioni numeriche
if ! [[ "$Threads" =~ ^[0-9]+$ ]] || (( Threads < 1 )); then
  echo "Errore: -t deve essere intero positivo." >&2; exit 2
fi

if ! [[ "$MaxIntronLength" =~ ^[0-9]+$ ]] || (( MaxIntronLength < 1 )); then
  echo "Errore: -m deve essere intero positivo." >&2; exit 2
fi

#####################################
# Env & software
#####################################

regtools_bin="/home/idris/software/regtools/build/regtools"

# proteggi gli hook con set -u

set +u
source /home/chiara/miniconda3/etc/profile.d/conda.sh
conda activate /home/chiara/miniconda3/envs/structural_annotation
set -u

#####################################
# Hisat2 index
#####################################

echo "[INFO] hisat2-build -> ${idx_Genome_Name}"
hisat2-build -p "$Threads" "$genome" "$idx_Genome_Name"

#####################################
# Hisat2 Single-End
#####################################

if [[ -n "$listsingle" ]]; then
  i=0
  while IFS= read -r r1; do
    [[ -z "$r1" || "$r1" =~ ^# ]] && continue
    ((i++))
    sample="$(basename "$r1")"; sample="${sample%%.*}"

    mkdir -p "RNAseq_SE_${i}"
    pushd "RNAseq_SE_${i}" >/dev/null
    echo "[SE] ${sample}"
   
hisat2 -p "$Threads" \
       -1 "$R1" \
       -2 "$R2" \
       -x "../${idx_Genome_Name}" \
       -k 1 \
| samtools view -b -@ "$Threads" > "${sample}.unsorted.bam"
       #--max-intronlen "$MaxIntronLength" \

    samtools sort -@ "$Threads" -o "${sample}.sorted.bam" "${sample}.unsorted.bam"
    samtools index -@ "$Threads" "${sample}.sorted.bam"
    rm -f "${sample}.unsorted.bam"
    popd >/dev/null
  done < "$listsingle"

fi

#####################################
# Hisat2 Paired-End (SEZIONE INVARIATA)
#####################################

if [[ -n "$listpaired" ]]; then
  i=0
  while IFS=$'\t' read -r R1 R2 extra; do
    [[ -z "${R1:-}" || "$R1" =~ ^# ]] && continue

    if [[ -z "${R2:-}" || -n "${extra:-}" ]]; then
      echo "Errore formato PE (serve: R1<TAB>R2). Riga: $R1 $R2 $extra" >&2
      exit 2
    fi

    [[ -f "$R1" ]] || { echo "File R1 non trovato: $R1" >&2; exit 2; }
    [[ -f "$R2" ]] || { echo "File R2 non trovato: $R2" >&2; exit 2; }

    ((i++))

    b="$(basename "$R1")"
    b="${b%.fastq.gz}"; b="${b%.fq.gz}"
    b="${b%.fastq}";   b="${b%.fq}"
    b="${b%_R1}"; b="${b%-R1}"; b="${b%.R1}"
    b="${b%_1}";  b="${b%-1}";  b="${b%.1}"
    sample="$b"

    echo "[INFO] hisat2 in paired end (PE) mode of ${b} sample"
   
    mkdir -p "RNAseq_PE_${i}"

    pushd "RNAseq_PE_${i}" >/dev/null
    echo "[PE] ${sample}"
   
    hisat2 -p "$Threads" \
           -1 "$R1" \
           -2 "$R2" \
           -x "../${idx_Genome_Name}" \
           -k 1 \
    | samtools view -b -@ "$Threads" > "${sample}.unsorted.bam"
          # --max-intronlen "$MaxIntronLength" \

    samtools sort -@ "$Threads" -o "${sample}.sorted.bam" "${sample}.unsorted.bam"

    samtools index -@ "$Threads" "${sample}.sorted.bam"

    rm -f "${sample}.unsorted.bam"

    popd >/dev/null

  done < "$listpaired"

fi


#####################################
# Merge & junctions -> AUGUSTUS hints
#####################################

# raccogli tutti i bam ordinati per coordinate

find . -maxdepth 2 -type f -path "./RNAseq_*/*.sorted.bam" | sort > bam.coord.list

if [[ ! -s bam.coord.list ]]; then
  echo "[ERRORE] Nessun BAM trovato in ./RNAseq_*/*.sorted.bam" >&2
  exit 2
fi

echo "[INFO] samtools merge"

samtools merge -@ "$Threads" -b bam.coord.list -o merged.coord.bam
samtools sort  -@ "$Threads" -o merged_sorted.bam merged.coord.bam
samtools index -@ "$Threads" merged_sorted.bam
echo "[INFO] regtools junctions extract"

#regtools
"$regtools_bin" junctions extract -o junctions.bed -s RF merged_sorted.bam


# ------------------------------------------------------------------

# Conversione junctions.bed (BED12) → GFF intron HINTS

# - Calcolo intron_start/end dai blockSizes (2 blocchi) di BED12:

#     intronStart0 = chromStart + blockSizes[0]

#     intronEnd0   = chromEnd   - blockSizes[last]

#   (BED è 0-based half-open; GFF è 1-based closed → start +1)

# - Attributi:

#     src=E           (Expression/RNA-seq)

#     grp=JUNC_chr_start_end_strand

#     mult=<score>    (#reads; colonna 5 BED)

# ------------------------------------------------------------------


awk -v OFS="\t" '

  BEGIN{ FS="\t" }

  /^#/ || NF<6 { next }  # serve almeno BED6
  {
    chr=$1; chromStart=$2; chromEnd=$3; name=$4; score=$5; strand=$6;

    # di norma regtools emette BED12 con due blocchi (ancore esoniche)
    # campi 10,11,12: blockCount, blockSizes, blockStarts

    blockCount = ($10=="" ? 0 : $10);
    blockSizes = $11;
   
    # fallback: se non è BED12 a 2 blocchi, salta (evita introni inconsistenti)

    if (blockCount != 2 || blockSizes == "") next;
    n=split(blockSizes, bs, ",");

    if (n<2 || bs[1]=="" || bs[2]=="") next;
    leftBlockLen  = bs[1];
    rightBlockLen = bs[2];

    intronStart0 = chromStart + leftBlockLen;
    intronEnd0   = chromEnd   - rightBlockLen;

    # controlli di coerenza

    if (intronEnd0 <= intronStart0) next;

    # converti a GFF (1-based closed): start = intronStart0+1, end = intronEnd0

    gffStart = intronStart0 + 1
    gffEnd   = intronEnd0;



    # costruisci grp stabile per questo introne

    grp = "JUNC_" chr "_" gffStart "_" gffEnd "_" strand;


    # se score vuoto, metti 1

    if (score == "" || score == ".") score=1;

    # stampa riga GFF intron

    if (score >= 10) print chr, "regtools", "intron", gffStart, gffEnd, ".", strand, ".", \
          "src=E;grp=" grp ";pri=4;mult=" score;
  }

' junctions.bed > RNAseq.unsorted.hints

# ordina per seqid,start,end,strand

sort -k1,1 -k4,4n -k5,5n -k7,7 RNAseq.unsorted.hints > HiSat2_min10.hints

# QA veloce

echo "[INFO] Tipi presenti in HiSat2_min10.hints:" >&2

awk -F"\t" '{print $3}' HiSat2_min10.hints | sort | uniq -c | sed -e "s/^/  /" >&2

# statistiche 
for bam in RNAseq_*/*.sorted.bam; do
    sample=$(basename "$bam" .sorted.bam)
    echo -n "$sample " 
    samtools flagstat "$bam" | awk 'NR==1 {total=$1} NR==2 {mapped=$1} END {printf "total:%d mapped:%d %.2f%%\n", total, mapped, mapped/total*100}'
done 

# pulizia
#rm -f merged.coord.bam RNAseq.unsorted.hints

echo "[DONE] Output: HiSat2_min10.hints"
