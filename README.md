# Lupinus_Genome_Annotation_Pipeline
This repository contains the modular and reproducible workflow used for the structural and functional annotation of Lupinus genomes. It use an Evidence-Driven Workflow for Structural Annotation Using Combined Data Sources from ab initio predictors, extrinsic biological evidence, and consensus modeling.

## Pipeline Architecture
The workflow is divided into two part, each containing a sequential list of scripts located in the scripts/ directory:

Sequential_Structural_Annotation.sh – From raw assembly and evidence to GFF3 consensus.

Functional_Annotation.sh – From protein FASTA to integrated functional TSV.


## Directory & Workflow Architecture
Lupinus_Annotation/
├── 1_Input_Genome/       # Raw and renamed assemblies (seqtk)
├── 2_RNAseq_Evidence/    # HISAT2 mappings, RegTools junctions, StringTie assemblies
├── 3_Protein_Evidence/   # Miniprot homology evidence
├── 4_Repeat_Masking/     # RepeatModeler library and RepeatMasker softmasked genome
├── 5_AbInitio_Models/    # BUSCO-trained Augustus configuration & Helixer inferences
├── 9_EVM/                # EvidenceModeler integration (weights.txt)
└── scripts/              # Executable bash scripts (0 to 10)


## Sequential Module Inventory
### 0. Directory Setup & Preprocessing
Scripts: Directory_Creation.sh / Directory_List.sh / rename_fasta.pl
Action: Generates standard directory architectures, logs structures, and strips complex chromosome headers using seqtk to avoid downstream pipeline crashes.

### 1. Repeat Identification and Genome Masking (Est. Run Time: 2–3 Days)
Script: 1_Repeat_Annotation.sh (calls /home/leonardo/annotation_script/repeat_to_gff.pl and convert_gff.sh)
Core Tools: RepeatModeler v2.0.6 $\rightarrow$ RepeatMasker v4.1.7p1.
Execution Logic: 
  1. RepeatModeler builds a de novo consensus library specific to the target assembly (consensi.fa.classified).
  2.RepeatMasker executes nucleotide-vs-consensus alignments using the rmblastn core engine.
Key Outputs:
  *.renamed_softmasked.fa (Primary genome input for Augustus execution, containing soft-masked lowercase intervals).
  *repeat.gff (Coordinates of genomic repeats)*
  *.temp.fa.tbl (Summary statistics table of repeat landscapes).

### 2. Species-Specific Predictive Training
Script: 2_BUSCO_training
Core Tool: BUSCO v5.8.3
Action: Captures single-copy orthologs from the specified lineage database, runs a native meta-training iteration for Augustus, and physically copies the learned statistical matrices into the active $AUGUSTUS_CONFIG_PATH/species/.

### 3. Transcript Expression Mapping (RNA-to-Genome)
Script: 3_HiSat2_custom.sh
Parameters enforced: --max-intronlen 100000 -k 1 (Maintains only the absolute best alignment per read).
Input Data: Non-directional Illumina PE libraries ($100\text{–}150\text{ bp}$, $Q>30$, $20\text{–}30\text{M}$ reads/library). Merges Giacomo's multi-tissue set (pod, root, leaf) and Leonardo's pangenome pools (leaf, genotypes G12873 & MIDAS).
Processing: samtools sort $\rightarrow$ samtools merge $\rightarrow$ regtools junctions extract -s RF.
Validation Threshold: Junctions are strictly filtered retaining only those with $>10$ uniquely mapping reads (Dong et al. 2023).
Key Output: HiSat2.hints (Formatted with mandatory src=E attributes for Augustus integration).

### 4. Transcriptome Assembly Generation
Script: 4_StringTie.sh
Core Tool: StringTie v3.0.0
Parameters enforced: Illumina default tracking mode (-L disabled), -m 200, -a 10, -j 1.
Output Conversion: Parses the default Stringtie.RAW.gtf through custom awk matrices, re-coding structural entries into exonpart tokens, and attaching ;src=E;grp=<transcript_id> annotations.

### 5. Cross-Species Proteomic Homology Mapping
Script: 5_Miniprot.sh
Core Tools: miniprot $\rightarrow$ miniprot-boundary-scorer $\rightarrow$ miniprothint.py
Reference inputs: Core high-quality proteomes (Glycine max, Phaseolus vulgaris, Lupinus angustifolius, Medicago truncatula).
Processing Workflow: Indexes target genome (genome.mpi) $\rightarrow$ Aligns proteins using dynamic programming heuristics $\rightarrow$ Filters structural shifts using a BLOSUM62.csv matrix filter.
Key Output: Proteome_LAST.hints (Standardized protein coordinate anchors embedded with src=P).

### 6. Intrinsic Gene Inferences
Script: 6_Augustus.sh
Action: Launches HMM-driven gene prediction models on the renamed_softmasked.fa genome. It dynamically ingests HiSat2.hints and Proteome_LAST.hints, reading relative validation coefficients through a custom extrinsic.cfg configuration file.

### 7. Structural Model Refining & External Inferences
Script: 7_interpro_filtering.sh (calls extract_interproscan.sh)
Action: Processes Augustus GFF models through InterProScan v5.72-103.0 to filter out computational fragments lacking real protein coding potential. In parallel, Helixer deep learning models are executed on the raw unmasked assembly.

### 8. Weighted Consensus Assembly
Script: 8_EVM_V6.sh
Core Tool: EvidenceModeler (EVM) v2.1.0
Action: Consolidates structural layers using an empirically sound weights.txt file configured for legume genome distributions.

### 9. Post-Consensus Diagnostics & Functional Deployment
Scripts: 9_Extract_data_From_GFF & 10_Functional_analysis_oneclick_V3.sh
Action: Segregates non-coding RNAs (ncRNAs), purges unexpressed single-exon detritus, and launches automated sequence homology lookups (DIAMOND BLASTp, EggNOG-mapper, and GO catalog allocations).
