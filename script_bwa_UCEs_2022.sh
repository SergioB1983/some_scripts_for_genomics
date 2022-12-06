########################
#### PASOS A SEGUIR ####
########################

# Until phyluce 1.6, there was a script that used the an unique reference to create bam files from the fasta data extracted after the 'phyluce_assembly_get_fastas_from_match_counts' and the 'phyluce_assembly_explode_get_fastas_file' scripts: 'phyluce_snp_bwa_multiple_align'

# The last version of phyluce (1.7) dont have this script, which is important to implementing the phasing processes (to create a SNPs matrix) as implemented by Harvey et al. (2016) and modified by Glaucia del Rio (Cornell Laboratory of Ornithology).

# With this in mind, I create a serie of commands that allowed to create the bam files necessaries to continue the phasing process of UCEs data.

# The codification in the bam files was following the script 'phyluce_snp_bwa_multiple_align' of the 1.6 version of phyluce.

## Preliminary step> create the dictionary for the sample that will be used as reference

# To create index of the reference sample with samtools
bwa index ./exploded-fastas/reference.unaligned.fasta
#
samtools faidx ./exploded-fastas/reference.unaligned.fasta

# To create .dict file of the reference (using the gatk software)
./gatk-4.3.0.0/./gatk CreateSequenceDictionary --REFERENCE ./exploded-fastas/reference.unaligned.fasta

## First step> Create the folders where you will do the phasing

# List the folders in the clean-read folder
cd ./clean-fastq/
#
find -type d -maxdepth 1 > ./clean-fastq/list_folders.txt

# Create and go to a new folder where you will do the phasing, and run this:
mkdir SNPs_phasing
#
cd ./SNPs_phasing
#
mkdir $( tr '[,\n]' ' ' < ./clean-fastq/list_folders.txt )

## Second step> Create the sam, bam files for each sample like the 'phyluce_snp_bwa_multiple_align' script

# Go to the SNPs folder created in the anterior step.
cd ./clean-fastq

# Using bwa-mem to create the intemedite sam file from the clean-reads with a unique reference (recursive)
# each sam file will be moved to its respective folder.
for dir in *; do
    for r in 1.fastq.gz; do
        for j in 2.fastq.gz; do
            outfile=${dir%/}.sam
            ref="./exploded-fastas/reference.unaligned.fasta"
            glob1=split-adapter-quality-trimmed/*${r}
            glob2=split-adapter-quality-trimmed/*${j}
            out_folder=./SNPs_phasing/$(basename -- "$dir")/
            bwa mem -t 4 $ref "$dir"/$glob1 "$dir"/$glob2 > "$dir/$outfile" &&
            mv "$dir"/$outfile $out_folder
        done
    done
done

# Now, go to the SNPs folder created in the anterior step.
cd ./SNPs_phasing

# Using samtools to transform the sam file in a bam file (recursive)
# samtools view!!!
# set the '@'' parameter with the available cores from your PC/cluster
for dir in */; do
    for r in .sam; do
        outfile=${dir%/}.bam
        glob=*${r}
        samtools -@4 view -bS "$dir"/$glob -o "$dir/$outfile"
    done
done

# Using samtools to 'sort' the sam file (recursive)
for dir in */; do
    for r in .bam; do
        outfile=${dir%/}.sorted.bam
        glob=*${r}
        samtools -@4 sort "$dir"/$glob -o "$dir/$outfile"
    done
done

# Using samtools to index the 'sorted.bam' file (recursive)
for dir in *; do
    for r in .sorted.bam; do
        glob=*${r}
        samtools -@4 index "$dir"/$glob
    done
done

# Using the 'CleanSam' script of the gatk software to clean the bam file (recursive)
for dir in *; do
    for r in .sorted.bam; do
        outfile=${dir%/}-CL.bam
        glob=*${r}
        ./gatk-4.3.0.0/./gatk CleanSam -I "$dir"/$glob -O "$dir/$outfile"
    done
done

# Extract the flowcell codes from the clean-reads (recursive)
# You need to run this command in the clean-reads folder (here the './clean-fastq/' folder)
# After this script, you have to manually edit the .flowcell file, and just left the flowcell code in the text.
for dir in *; do
    for r in -READ1.fastq.gz; do
        outfile=${dir%/}.flowcell
        glob=split-adapter-quality-trimmed/*${r}
        out_folder=./SNPs_phasing/$(basename -- "$dir")/
        zcat "$dir"/$glob | sed -n '1p' > "$dir/$outfile" &&
        mv "$dir"/$outfile $out_folder
    done
done

# Using the 'AddOrReplaceReadGroups' script of the gatk software on the bam file (recursive)
for dir in *; do
    for r in -CL.bam; do
        outfile=${dir%/}-CL-RG.bam
        glob=*${r}
        flowcell=$(find $dir -type f -name '*.flowcell')
        code_flowcell=$(head -n1 $flowcell)
        ./gatk-4.3.0.0/./gatk AddOrReplaceReadGroups -I "$dir"/$glob -O "$dir/$outfile" --RGID "$dir" --RGLB 'Lib1' --RGPL 'illumina' --RGPU $code_flowcell --RGSM "$dir"
    done
done

# Using the 'MarkDuplicatesSpark' script of the gatk software on the bam file to mark duplicates (recursive)
# the option --spark-master local will be set the number of cores used.
for dir in *; do
    for r in -CL-RG.bam; do
        glob=*${r}
        outfile=${dir%/}-CL-RG-MD.bam
        stats=$(basename -- "$dir").txt
        ./gatk-4.3.0.0/./gatk MarkDuplicatesSpark -I "$dir"/$glob -O "$dir/$outfile" -M "$dir/$stats" --spark-master local[4]
    done
done

# Using the 'SetNmMdAndUqTags' script of the gatk software to fixes the NM, MD, and UQ tags in a SAM/BAM/CRAM file (recursive)
# this command is based in the Faircloth-lab documentation (https://protocols.faircloth-lab.org/_/downloads/en/latest/pdf/)
for dir in *; do
    for r in -RG-MD.bam; do
        glob=*${r}
        outfile=${dir%/}-CL-RG-MD.bam
        ref="./exploded-fastas/reference.unaligned.fasta"
        ./gatk-4.3.0.0/./gatk SetNmMdAndUqTags -I "$dir"/$glob -O "$dir/$outfile" --REFERENCE_SEQUENCE $ref
    done
done

# USing the 'BuildBamIndex' script to generates a BAM index ".bai" file.
for dir in *; do
    for r in -MD.bam; do
        glob=*${r}
        outfile=${dir%/}-CL-RG-MD-fix.bai
        ./gatk-4.3.0.0/./gatk BuildBamIndex -I "$dir"/$glob -O "$dir/$outfile"
    done
done

## From here, you can delete all intermetiate files, and use the '-CL-RG-MD.bam' files (as the 1.6 protocol of phyluce).

## The next step will be using the 'MergeSamFiles' script to merge the -MD.bam files and continues with the protocol of phasing.

## Bibliography

Faircloth BC. 2015. PHYLUCE is a software package for the analysis of conserved genomic loci. Bioinformatics. doi: 10.1093/bioinformatics/btv646.

Harvey, M. G., Smith, B. T., Glenn, T. C., Faircloth, B. C., Brumfield, R. T. (2016). Sequence Capture versus Restriction Site Associated DNA Sequencing for Shallow Systematics. Systematic Biology, 65(5), 910–924. https://doi.org/10.1093/sysbio/syw036.

