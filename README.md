# 2025_HRR_methods_CCR_HS

## Calling somatic variants

The following code is used to run GATK Mutect2.

```bash
# Running GATK Mutect2
$GATK --java-options "-Xmx6g" Mutect2 \
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
$GATK --java-options "-Xmx6g" LearnReadOrientationModel \
	 -I ${path_to_output_directory}/${sample_pair_name}.f1r2.tar.gz \
	 -O ${path_to_output_directory}/${sample_pair_name}.read-orientation-model.tar.gz

# Contamination calculation and filtering
$GATK GetPileupSummaries \
    -R ${path_to_fasta_file} \
    -I ${path_to_tumor_bam} \
    -O ${path_to_output_directory}/${tumor_sample_name}_getpileupsummaries.table

$GATK GetPileupSummaries \
    -R ${path_to_fasta_file} \
    -I ${path_to_normal_bam} \
    -O ${path_to_output_directory}/${normal_sample_name}_getpileupsummaries.table

$GATK CalculateContamination \
    -I ${path_to_output_directory}/${tumor_sample_name}_getpileupsummaries.table \
    -matched ${path_to_output_directory}/${normal_sample_name}_getpileupsummaries.table \
    -O ${path_to_output_directory}/${sample_pair_name}_calculatecontamination.table \
    -tumor-segmentation ${path_to_output_directory}/${sample_pair_name}_segments.table

$GATK FilterMutectCalls \
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
# Running Manta
python2 ${path_to_Manta}/configManta.py \
--normalBam ${path_to_normal_bam} \
--tumorBam ${path_to_tumor_bam} \ 
--referenceFasta ${path_to_fasta_file} \
--runDir ${path_to_output_directory} \

# Running Strelka2
python2 ${path_to_Strelka}/configureStrelkaSomaticWorkflow.py \
--normalBam ${path_to_normal_bam} \
--tumorBam ${path_to_tumor_bam} \
--ref ${path_to_fasta_file} \
--indelCandidates ${path_to_Manta_output_directory}/candidateSmallIndels.vcf.gz \
--runDir ${path_to_output_directory}

```

## Calling germline variants

The following code is used to run GATK HaplotypeCaller.

```bash
gatk HaplotypeCaller \
   -R ${path_to_fasta_file} \
   -I ${path_to_normal_bam} \
   -O ${path_to_output_directory}/${sample_name}.g.vcf.gz \
   -ERC GVCF \
   -L ${path_to_bed_file}
```

## Calling copy number variations

The following code is used to run FACETS.
```bash
# 1. Script----------
# script name = FACETS.py

# Read sample information
try:
    fp = open(sys.argv[argn])
    for line in fp.readlines():
        line = line.strip()
        if len(line) == 0:
            continue

        fields = line.split(" ")
        if len(fields) != 3:
            continue

        samples.append(fields[0])

        type = fields[2].lower()
        if type != 't' and type != 'n':
            continue

        sys.stdout.write("%s %s %s\n" % (fields[0], fields[1], type))

        if fields[1] not in pairs:
            pairs[fields[1]] = type
            pairsamples[fields[1]] = [fields[0]]
        else:
            pairs[fields[1]] += type
            pairsamples[fields[1]].append(fields[0])

    fp.close()
except IOError:
    pass

sys.stdout.write("EOF\n")

# Run FACETS
module load HTSlib/1.21-GCC-12.2.0
${path to FACETS}/inst/extcode/snp-pileup \
-g -q15 -Q20 -P100 -r25,0 \
${path to FACETS}/inst/extdata/common_all_20180418.vcf.gz \
${path_to_output_directory}/${sample_pair_name}.csv.gz \
${path_to_normal_bam} \
${path_to_tumor_bam}

module load R/4.2.0-foss-2020b
Rscript ${path to FACETS}/bin/cnv_facets.R -p ${path_to_output_directory}/${sample_pair_name}.csv.gz -o "${path_to_output_directory}/${sample_pair_name}_cnv_facets"

# 2. Example of input file----------
# input file name = samples.txt
[tumor_sample_dir_name] [sample_pair_name] [T(tumor)/N(normal)]
[normal_sample_dir_name] [sample_pair_name] [T(tumor)/N(normal)]

# 3. Running FACETS----------
python FACETS.py
```


## Calculating Genomic Instability Score using scarHRD
The following code is used to generate scarHRD input file.
```bash
# 1. Script----------
# script name = scarHRD_input.R

#!/usr/bin/Rscript
# Parse command line arguments
args = commandArgs(TRUE)

# Assign arguments to variables
subset = args[1]
loc = args[2]

# Load required libraries
library("facets")
library("data.table")

# Load FACETS output file
datafile = file.path(${path_to_FACETS_output_directory})

# Process the sample data
rcmat = readSnpMatrix(datafile)
xx = preProcSample(rcmat)
oo = procSample(xx, cval=150)
fit = emcncf(oo)

# Prepare scarHRD input data
seg <- data.frame(SampleID = subset,
                 Chromosome = fit$cncf$chrom,
                 Start_position = fit$cncf$start,
                 End_position = fit$cncf$end,
                 total_cn = fit$cncf$tcn.em,
                 A_cn = fit$cncf$tcn.em - fit$cncf$lcn.em,
                 B_cn = fit$cncf$lcn.em,
                 ploidy = fit$ploidy)

# Write scarHRD input file
write.table(seg,
            file = output_file,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)
            
# 2. Running script----------
Rscript scarHRD_input.R ${sample_pair_name} ${path_to_FACETS_output_directory}
```

The following code is used to run scarHRD.

```bash
# 1. Script----------
# script name = scarHRD.R

#!/usr/bin/Rscript
# Parse command line arguments
args = commandArgs(TRUE)

# Assign arguments to variables
subset = args[1]
loc = args[2]

# Load required library
library(scarHRD)

# Process the data and handle potential NA values
input_data = read.table(${path_to_scarHRD_input_file}, header=TRUE)

# Calculate scarhrd scores
tryCatch({
    scores = scar_score(input_data, 
                       reference = "grch38", 
                       seqz = FALSE)
    
    # Prepare output data
    result = data.frame(
        Sample = subset,
        HRD = scores[1],
        Telomeric_AI = scores[2],
        LST = scores[3],
        HRD_sum = scores[4]
    )
    
    # Write results to file
    write.table(result, 
                file = output_file, 
                quote = FALSE, 
                sep = "\t", 
                row.names = FALSE)

# 2. Running scarHRD----------
Rscript scarHRD.R ${sample_pair_name} ${path_to_scarHRD_input_directory}

```

## Identifying Mutational Signatures

The following code is used to run SigProfilerAssignment.

```bash
# Analysis environment setup----------
conda create --name sigprofiler -y
conda activate sigprofiler
conda install python=3.10 r-base r-devtools r-reticulate -c conda-forge -y
pip install SigProfilerMatrixGenerator numpy==1.23.5
pip install SigProfilerAssignment
pip install SigProfilerPlotting

# Running SigProfilerAssignment----------
conda activate sigprofiler
python
Analyze.cosmic_fit(samples=${path_to_vcf_directory}, output=${path_to_output_directory}, input_type="vcf", context_type="96", genome_build="GRCh38", cosmic_version=3.4)
```

## Reference
