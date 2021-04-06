STRAIN=Mt_h37rv

# Create the first vcf (all default parameters, brut-force writting of the parameters in the script)
./createIndividuals.py Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa  Mt_h37rv.vcf > run.val
# Annotation and removal of the creation of stop codon (I put it here in case we wanted to use snpEff that would have pull these one appart)
bcftools csq -f Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa -g Mycobacterium_tuberculosis_h37rv.ASM19595v2.46.chromosome.Chromosome.gff3 Mt_h37rv.vcf  > Mt_h37rv.annot.vcf 
grep -v stop  Mt_h37rv.annot.vcf >  Mt_h37rv.corr.vcf

picard CreateSequenceDictionary R=Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa O=Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.dict

# Add ##INFO=<ID=BCSQ,Number=1,Type=String,Description="BCSQ"> to  Mt_h37rv.corr.vcf
gatk IndexFeatureFile -I Mt_h37rv.corr.vcf


# And here it is: the creation of fastqs
I=`bcftools query -l ${STRAIN}.corr.vcf`
for sample in $I
do
    #sample=00RifRes121211201
    # Select the sample
    bcftools view -s ${sample} ${STRAIN}.corr.vcf > ${STRAIN}.${sample}.vcf &&
    
    # Remove the non-variant (Maybe optional)
    gatk SelectVariants --exclude-non-variants --variant ${STRAIN}.${sample}.vcf --output ${STRAIN}.${sample}.varOnly.vcf &&
    
    # Create the Fasta
    gatk FastaAlternateReferenceMaker -R=Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa -O ${sample}.fa --variant ${STRAIN}.${sample}.varOnly.vcf &&
    
    # Create the Fastqs (*.mut should be empty)
    wgsim -N 50000 -1 150 -2 150 -r 0 -R 0 ${sample}.fa ${sample}.forward.fq ${sample}.reverse.fq > ${sample}.mut &&
    rm ${STRAIN}.${sample}.varOnly.vcf*  ${sample}.fa ${STRAIN}.${sample}.vcf
done


bcftools +fill-tags ${STRAIN}.corr.vcf > ${STRAIN}.corr.annot.vcf
