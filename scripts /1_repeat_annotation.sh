#!/usr/bin/env bash

#####################################

# Usage & parsing

#####################################

usage() {

  cat <<'EOF'

Repeat annotation 28_8_2025 (f.fontana)

Questo script individua identifica le ripetizioni in un genoma e ne effettua il mascheramento.

USO:
  ./comando.sh -x <path/to/genome> -t threads -p project_name


Opzioni:

        -x (obbligatorio): path al genoma da analizzare e mascherare;
        -t (opzionale): numero di threads da usare (default 50).
        -p (opzionale): nome del progetto (default Project_One)

EOF

}

  

genome=""
nume_threads=""
project_name=""

while getopts ":x:t:p:h" opt; do
  case "${opt}" in
    x) genome="${OPTARG}" ;;
    t) nume_threads="${OPTARG}" ;;
    p) project_name="${OPTARG}" ;;
    h) usage; exit 0 ;;
    \?) echo "Errore: opzione non valida -${OPTARG}" >&2; usage; exit 2 ;;
    :)  echo "Errore: l'opzione -${OPTARG} richiede un argomento." >&2; usage; exit 2 ;;
  esac
done

shift $((OPTIND-1))

if [[ -z "${genome}" ]]; then
  echo "Errore: Il path genoma (-x) è obbligatorio." >&2
  usage
  exit 2
fi

if ! [[ "${nume_threads}" =~ ^[0-9]+$ ]]; then
  echo "Errore: -t deve essere un numero intero positivo" >&2
  exit 2
fi

if [[ -z "${nume_threads}" ]]; then
  echo "Threads non specificati. Verranno usati 50 threads as default" >&2
  nume_threads=50
fi

if ! [[ "${nume_threads}" =~ ^[0-9]+$ ]]; then
  echo "Errore: -t deve essere un numero intero positivo" >&2
  exit 2
fi

if [[ -z "${project_name}" ]]; then
  echo "Nome del progetto non specificato. Verranno usato <Project_One> as default" >&2
  project_name="Project_One"
fi

  

#####################################

# Main block

#####################################

  

source /home/chiara/miniconda3/etc/profile.d/conda.sh
conda activate /home/chiara/miniconda3/envs/repeatmodeler

# Stop the script if any command fails
set -e

seqtk rename $genome seq > temp.fa 
# perdi definitivamente i nomi originali !!! REMEMBER 

# tmp.fa genoma rinominato in <seq1 ...

mkdir Database_${project_name};
cd Database_${project_name};
BuildDatabase -name $project_name ../temp.fa
cd ..;

RepeatModeler -threads $nume_threads -database Database_${project_name}/${project_name} -LTRStruct temp.fa

# Better way to find and copy the library
CONSENSUS=$(find RM_* -name "consensi.fa.classified" | head -n 1)
if [[ -z "$CONSENSUS" ]]; then
    echo "ERRORE: consensi.fa.classified non trovato"
    exit 1
fi
cp "$CONSENSUS" ./consensi.fa.classified

# controllo esistenza consensi.fa
if [[ ! -s consensi.fa.classified ]]; then
  echo "ERRORE: consensi.fa.classified non trovato o vuoto"
  exit 1
fi

## si è interrotto qui perchè non trovava consensi.fa visto che non era riuscito a copiarlo 

RepeatMasker -pa $nume_threads -xsmall -lib ./consensi.fa.classified temp.fa

# Se usa HMM forzare ncbi
#RepeatMasker - engine ncbi -pa $nume_threads -xsmall -lib ./consensi.fa.classified temp.fa

# Uso Dfam_5 Viridipiante come libreria 
#RepeatMasker -pa $nume_threads -xsmall -species Viridiplantae temp.fa

# -xsall returns repetitive regions in lowercase rather than masked (softmask)

mv temp.fa.masked genome_masked.fa

# L'intestazione originale viene riportata in repeat.gff

perl /home/chiara/Scripts/Annotation/repeat_to_gff.pl temp.fa.out

grep '>' temp.fa | sed 's/>//g' > list1
grep '>' $genome | sed 's/>//g' > list2
paste list1 list2 -d , > list

/home/chiara/Scripts/Annotation/convert_gff.sh list temp.fa.out.gff > repeat.gff

rm -f list
rm -f list1
rm -f list2
  
# L'intestazione originale di $masked_genome viene sostituita con le intestazioni di $genome e fixato.

less genome_masked.fa | grep ">" | sed 's/>//g' > a
less ${genome} | grep ">" | sed 's/>//g' > b
paste -d "\t" a b > list

perl /home/chiara/Scripts/Annotation/rename_fasta.pl list genome_masked.fa > renamed_softmasked.fa

sed -i -r 's/^(>[^ ]+).*/\1/' renamed_softmasked.fa

rm -f a
rm -f b
rm -f list
