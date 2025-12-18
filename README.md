# 2025_HRR_methods_CCR_HS

Tumor samples were classified as homologous recombination (HR)-deficient if they met at least one of the following criteria:

1. Genomic Instability Score (GIS) ≥ 42
2. Presence of a pathogenic or likely pathogenic somatic or germline mutation in BRCA1 or BRCA2

Samples that did not meet either of these criteria were classified as HR-proficient.

## Variant calling
### 1. Detection of somatic SNVs and INDELs
Somatic SNVs were detected using GATK Mutect2, and somatic INDELs were identified based on the overlap between Mutect2 and Strelka2.

The following code is used to run GATK Mutect2.

```bash
# GATK version : 4.6.0.0

# Running GATK Mutect2
# Example of fasta file : hs38DH.fasta
# Example of a panel of normals (PON) file: 1000g_pon.hg38.vcf.gz

gatk Mutect2 \
        -R ${path_to_fasta_file} \
        -I ${path_to_tumor_bam} \
        -I ${path_to_normal_bam} \
        -tumor ${tumor_sample_name} \
        -normal ${normal_sample_name} \
        --panel-of-normals ${path_to_PON_file} \
        --f1r2-tar-gz ${path_to_output_directory}/${sample_pair_name}.f1r2.tar.gz \
        -O ${path_to_output_directory}/${sample_pair_name}.mutect2.vcf \
        -bamout ${path_to_output_directory}/${sample_pair_name}.mutect2.bam

# Running GATK LearnReadOrientationModel
gatk LearnReadOrientationModel \
	 -I ${path_to_output_directory}/${sample_pair_name}.f1r2.tar.gz \
	 -O ${path_to_output_directory}/${sample_pair_name}.read-orientation-model.tar.gz

# Contamination calculation and filtering
gatk GetPileupSummaries \
    -R ${path_to_fasta_file} \
    -I ${path_to_tumor_bam} \
    -O ${path_to_output_directory}/${tumor_sample_name}_getpileupsummaries.table

gatk GetPileupSummaries \
    -R ${path_to_fasta_file} \
    -I ${path_to_normal_bam} \
    -O ${path_to_output_directory}/${normal_sample_name}_getpileupsummaries.table

gatk CalculateContamination \
    -I ${path_to_output_directory}/${tumor_sample_name}_getpileupsummaries.table \
    -matched ${path_to_output_directory}/${normal_sample_name}_getpileupsummaries.table \
    -O ${path_to_output_directory}/${sample_pair_name}_calculatecontamination.table \
    -tumor-segmentation ${path_to_output_directory}/${sample_pair_name}_segments.table

gatk FilterMutectCalls \
    -V ${path_to_output_directory}/${sample_pair_name}.mutect2.vcf.gz \
    -R ${path_to_fasta_file} \
    --ob-priors ${path_to_output_directory}/${sample_pair_name}.read-orientation-model.tar.gz \
    --contamination-table ${path_to_output_directory}/${sample_pair_name}_calculatecontamination.table \
    --tumor-segmentation ${path_to_output_directory}/${sample_pair_name}_segments.table \
    --stats ${path_to_output_directory}/${sample_pair_name}.mutect2.vcf.gz.stats \
    -O ${path_to_output_directory}/${sample_pair_name}.mutect2.filtered.vcf
```

The following code is used to run Manta and Strelka2.

```bash
# Manta version : 1.6.0
# Sterlka2 version : 2.9.10
# configManta.py and configureStrelkaSomaticWorkflow.py are available upon installation of Manta and Strelka.

# Running Manta
python2 manta/configManta.py \
--normalBam ${path_to_normal_bam} \
--tumorBam ${path_to_tumor_bam} \ 
--referenceFasta ${path_to_fasta_file} \
--runDir ${path_to_output_directory} \

# Running Strelka2
python2 strelka/configureStrelkaSomaticWorkflow.py \
--normalBam ${path_to_normal_bam} \
--tumorBam ${path_to_tumor_bam} \
--ref ${path_to_fasta_file} \
--indelCandidates ${path_to_Manta_output_directory}/candidateSmallIndels.vcf.gz \
--runDir ${path_to_output_directory}
```

### 2. Detection of germline SNVs and INDELs

The following code is used to run GATK HaplotypeCaller.

```bash
gatk HaplotypeCaller \
   -R ${path_to_fasta_file} \
   -I ${path_to_normal_bam} \
   -O ${path_to_output_directory}/${sample_name}.g.vcf.gz \
   -ERC GVCF
```

### 3. Detection of somatic CNVs
The following code is used to run FACETS.
```bash
# FACETS version : 0.6.2
# snp-pileup and cnv_facets.R are available upon installation of the FACETS.

facets/inst/extcode/snp-pileup \
${path_to_output_directory}/${sample_pair_name}.csv.gz \
${path_to_normal_bam} \
${path_to_tumor_bam}

Rscript facets/bin/cnv_facets.R -p ${path_to_output_directory}/${sample_pair_name}.csv.gz -o ${path_to_output_directory}/${output_file}
```

## Downstream analyses for HRD determination

### 1. Calculating GIS
The following code is used to run scarHRD.
```bash
# R environment
# scarHRD version : 0.1.1

library(scarHRD)
scores <- scarHRD::scar_score(input_data)
```

### 2. Identifying Mutational Signatures
The following code is used to run SigProfilerAssignment.
```bash
# SigProfilerAssignment version : 0.1.9

pip install SigProfilerAssignment
pip install SigProfilerMatrixGenerator
```
Execute the following code in Python.
```bash
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh38')
from SigProfilerAssignment import Analyzer as Analyze
Analyze.cosmic_fit(${path_to_vcf_directory}, ${path_to_output_directory}, input_type="vcf", context_type="96", genome_build="GRCh38", cosmic_version=3.4)
```

## Reference
1. McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010;20(9):1297–1303.
2. Kim S, Scheffler K, Halpern AL, et al. Strelka2: fast and accurate calling of germline and somatic variants. Nat Methods. 2018;15(8):591–594.
3. Shen R, Seshan VE. FACETS: allele-specific copy number and clonal heterogeneity analysis tool for high-throughput DNA sequencing. Nucleic Acids Res. 2016;44(16):e131.
4. Sztupinszki Z, Diossy M, Krzystanek M, et al. Migrating the SNP array-based homologous recombination deficiency measures to next generation sequencing data of breast cancer. NPJ Breast Cancer. 2018;4:16.
5. Díaz-Gay M, Vangara R, Barnes M, et al. Assigning mutational signatures to individual samples and individual somatic mutations with SigProfilerAssignment. Bioinformatics. 2023;39(12).



