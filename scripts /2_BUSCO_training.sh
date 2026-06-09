#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Busco training of Augustus 28_8_2025
USO:  ./comando.sh -x <genome.fa> -t <threads> -s <species> -l <lineage>
EOF
}

# default threads
Threads=50
nume_threads=""
species=""
genome=""
Lineage=""

while getopts ":x:t:s:l:h" opt; do
  case "${opt}" in
    x) genome="${OPTARG}" ;;
    t) nume_threads="${OPTARG}" ;;
    s) species="${OPTARG}" ;;
    l) Lineage="${OPTARG}" ;;
    h) usage; exit 0 ;;
    \?) echo "Errore: opzione non valida -${OPTARG}" >&2; usage; exit 2 ;;
    :) echo "Errore: l'opzione -${OPTARG} richiede un argomento." >&2; usage; exit 2 ;;
  esac
done
shift $((OPTIND-1))

# validazioni
[[ -s "$genome" ]] || { echo "Errore: genoma non trovato" >&2; exit 2; }
[[ -n "$species" ]] || { echo "Errore: nome specie mancante" >&2; exit 2; }
[[ -n "$Lineage" ]] || { echo "Errore: lineage BUSCO mancante" >&2; exit 2; }

# threads
if [[ -n "$nume_threads" ]]; then
  if ! [[ "$nume_threads" =~ ^[0-9]+$ ]]; then
    echo "Errore: -t deve essere numero intero positivo" >&2; exit 2
  fi
  Threads=$nume_threads
fi
echo "[INFO] USO $Threads threads"

#####################################
# Conda + Java memory
#####################################
source /home/chiara/miniconda3/etc/profile.d/conda.sh
conda activate /home/chiara/miniconda3/envs/busco 

# export della heap per BBTools stats.sh
export BBMAP_JAVA_OPTS="-Xmx64g"

#####################################
# Esecuzione BUSCO
#####################################

busco -i "$genome" -l "$Lineage" -m genome -o "$species" -c "$Threads" --augustus 


# copia i parametri per Augustus
cp -r "$species/run_$Lineage/augustus_output/retraining_parameters/BUSCO_$species" \
   /home/chiara/miniconda3/envs/Augustus/config/species/
