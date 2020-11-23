### CRCmapper: MAP CORE REGULATORY CIRCUITRY 

REFERENCE:
Models of Human Core Transcriptional Regulatory Circuitries. Violaine Saint-André, Alexander J. Federation, Charles Y. Lin, Brian J. Abraham, Jessica Reddy, Tong Ihn Lee, James E. Bradner, Richard A. Young. Genome Res. 2016. 26: 385-396 (<https://genome.cshlp.org/content/26/3/385.long>)

Please cite this article when using this code.

SOFTWARE AUTHORS: Violaine Saint-Andre, Alexander J. Federation, Charles Y. Lin.

CONTACT: <violaine.saint-andre@curie.fr>

Developed using `Python 2.7.3`

PURPOSE:
To build Core Regulatory Circuitry from H3K27ac ChIP-seq data

#### 1. REQUIREMENTS
FIMO (Grant et al. 2011) from the MEME suite and SAMtools (http://www.htslib.org/) (Li et al., 2009) must be installed

Code must be run from the directory in which it is stored together with `TFlist_NMid_hg.txt`, `TFlist_NMid_ms.txt`, `VertebratePWMs.txt`, `MotifDictionary.txt`, `bamToGFF.py`, `utils.py`, `bamToGFFutils.py` and a directory named "annotation" containing the genome annotation files (`hg19_refseq.ucsc`, `hg18_refseq.ucsc`, `mm9_refseq.ucsc`)

The bam file of sequencing reads for H3K27ac must be sorted and indexed using SAMtools

Fasta files for the genome used must be placed in a directory that will be specified when runing the program (`-f` option). Those files must be split by chromosome and termed "chrN.fa" with N being the chromosome number. They can be downloaded from <http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/> (they will need to be unzipped)

#### 2. CONTENT

`CRCmapper.py`: main program

`utils.py`: utility methods

`TFlist_NMid_hg.txt`: TFs used and their human NMIDs

`TFlist_NMid_ms.txt`: TFs used and their murine NMIDs

`VertebratePWMs.txt`: Vertebrate Motifs PWM library

`MotifDictionary.txt`: TFs used and their associated motif PWM names

`bamToGFF.py`: program calculating density of sequencing reads from the bam file in specified regions, and the genome annotation file (<https://github.com/bradnerComputation/pipeline/blob/master/bamToGFF.py>) (Lin et al. 2012)

`annotation/hg19_refseq.ucsc`, `annotation/hg18_refseq.ucsc` and `annotation/mm9_refseq.ucsc`: genome annotation files

#### 3. USAGE

The program is run by calling CRCmapper.py from the directory containing all the documents:

`python CRCmapper.py -e [ENHANCER_FILE] -b [BAM_FILE] -g [GENOME] -f [FASTA]` 
`-s [SUBPEAKS] -x [EXP_CUTOFF] -l [EXTENSION-LENGTH] -n [NAME] -o [OUTPUT_FOLDER] [optional: ]`

Required parameters:

`-e` [ENHANCER_FILE]

[ENHANCER_FILE]: enhancer table (SAMPLE_AllEnhancers.table.txt) generated with ROSE (https://bitbucket.org/young_computation/rose)

`-b` [BAM_FILE]

[BAM_FILE]: sorted indexed bam file for H3K27ac sequencing reads

`-g` [GENOME]

[GENOME]: build of the genome to be used for the analysis. Currently supports HG19, HG18 and MM9

`-f` [FASTA]

[FASTA]: path to the corresponding genome fasta files

`-s` [SUBPEAKS]

[SUBPEAKS]: bedfile of peaks output from MACS used to identify SE constituents

`-x` [EXP_CUTOFF]

[EXP_CUTOFF]: percentage of transcripts that are not considered expressed, default=33

`-l` [EXTENSION-LENGTH]

[EXTENSION-LENGTH]: length (in bp) to extend constituents for motif search, default=500

`-n` [NAME]

[NAME]: sample name

`-o` [OUTPUT_FOLDER]

[OUTPUT_FOLDER]: directory to be used for storing output

Optional parameters:

`-a` [ACTIVITY_TABLE]

[ACTIVITY_TABLE]: a two column table with refseq in the first column and the associated activity (expression or promoter acetylation level) in the second column

`-E` [ENHANCER_NUMBER]

[ENHANCER_NUMBER]: the number of top ranked enhancers to include in the analaysis, default = supers

#### 4. OUTPUT FILES

`fimo.txt`: output of the motif search from the FIMO program

`matrix.gff`: location matrix in gff format used by bamToGFF.py program

`SAMPLE_ASSIGNMENT_GENES.txt`: list of gene names for genes assigned to SEs

`SAMPLE_ASSIGNMENT_TRANSCRIPTS.txt`: Transcripts NMIDs for transcripts assigned to SEs

`SAMPLE_AUTOREG.txt`: list of TFs gene names predicted to bind their own SE

`SAMPLE_bg.meme`: DNA background sequence file used with FIMO

`SAMPLE_CANDIDATE_TF_AND_SUPER_TABLE.txt`: table containing the candidate TFs and the location of their associated SEs

`SAMPLE_CRC_SCORES.txt`: all possible CRCs, ranked based on the average frequency of occurrence of the TFs they contain across all the possible interconnected auto regulatory loops (see reference above for details)

`SAMPLE_EXPRESSED_GENES.txt`: list of genes considered expressed (see reference above for details)

`SAMPLE_EXPRESSED_TRANSCRIPTS.txt`: list of transcripts considered expressed as explained in the reference above

`SAMPLE_subpeaks.bed`: bedfile of SE constituent sequences

`SAMPLE_SUBPEAKS.fa`: fasta file of SE constituent sequences used with FIMO

`SAMPLE_TSS.gff`: gff file with TSS coordinates used by bamToGFF.py

`TF_SAMPLE_motifs.bed`: DNA binding motif locations in extended enhancer constituents for this TF

**References**

Grant CE, Bailey TL, Noble WS. 2011. FIMO: scanning for occurrences of a given motif. Bioinformatics 27: 1017–8.

Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R. 2009. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25: 2078–9.

Lin CY, Lovén J, Rahl PB, Paranal RM, Burge CB, Bradner JE, Lee TI, Young RA. 2012. Transcriptional Amplification in Tumor Cells with Elevated c-Myc. Cell 151: 56–67.
