
### Table of contents

<details><summary>See</summary>  


<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:0 orderedList:0 -->

- [Table of contents](#table-of-contents)
- [Overview](#overview)
- [Background](#background)
- [Day1](#day1)
	- [Objectives](#objectives)
	- [Task 1:](#task-1)
  - [data day 1](#data-day-1)
	- [Task 2:](#task-2)
	- [Task 3:](#task-3)
	- [Task 4:](#task-4)
	- [Wrap up day 1](#wrap-up-day-1)
- [Day 2](#day-2)
	- [Objectives](#objectives)
	- [Data day 2](#data-day-2)
	- [Task 1](#task-1)
	- [Task2](#task2)
	- [Task3](#task3)
	- [Task4](#task4)
	- [Task5](#task5)
	- [Discussion and Wrap-up.](#discussion-and-wrap-up)

<!-- /TOC -->
</details>  


### Overview
![visual overview](figures/UU-NGS02_v2.svg)

### Background
Today and tomorrow we will work with **simulated** sequencing data derived from [_Mycobacterium tuberculosis_](https://en.wikipedia.org/wiki/Mycobacterium_tuberculosis).
 _Mycobacterium tuberculosis_ is a pathogenic bacteria, and the causative agent of [tuberculosis](https://en.wikipedia.org/wiki/Tuberculosis), a bacterial disease that kills about [1.2 million people annually](https://www.who.int/tb/publications/global_report/en/), making it one of the leading cause of death among infectious diseases.
 Since the treatment of tuberculosis requires sustained use of antibiotics, the rise of  multiple drug-resistant tuberculosis (MDR-TB) and extensively drug-resistant tuberculosis (XDR-TB) provides a serious concern for world health.

We will work with several "patient derived samples": strains that have been either sampled from patients where the disease responded well to antibiotic treatment, or from patients with a more problematic disease progression, where the infection did not respond to first-line antibiotics.

## Day1  

### Objectives
 - practise/re-use what you've learned yesterday on a different dataset. If anything was unclear yesterday, now is the time to catch up on that.
 - get a vague insight into how different data requires different treatment.
 - transition [_"from copy-pasting magic incantations"_ to _"using tools that i can adapt to my data"_.](modular.md)
 - gain increased familiarity with some of the basic data formats and handling them, using simple command line tools.

### Task 1:  
SSH into Rackham and request a interactive session, like you did yesterday.

<details><summary>tips</summary>

<p>

```bash
ssh -Y <user_name>@rackham.uppmax.uu.se
interactive -A UPPMAX 2023-2-24 -M snowy -t 4:00:00 --reservation=uppmax2023-2-24_7
```
</p>
</details>


<br>
<br>

## data day 1
We have prepared a directory with **.fastq** files from **100** samples:  
```/proj/g2020004/private/computer_practicals/NGS_workflow_day3_4/data ```

each sample consists of two files, containing the forward and reverse reads.
They follow the naming scheme   ```SampleID.DIRECTION.FILEENDING```.
For example, sample **01200211Eth2Res120** , will look like this:

```
01200211Eth2Res120.forward.fq
01200211Eth2Res120.reverse.fq
```
In the same folder, you will also find the reference file:

```Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna.chromosome.Chromosome.fa```
which you will need for alignment.


### Task 2:
pick two samples (so, 4 files) and [copy](https://linux.die.net/man/1/cp) or [simlink](https://ss64.com/bash/ln.html) them to your home directory.

**Questions:**

 - _How does this data differ from the one that you worked on before?_
 -  _Do you think you will be able to use the GATK BSQR tool on this dataset?_
 - _Can you figure out how many reads there are in one file?_
<details><summary>tips</summary>
    <p>

 [ the fastq format](https://en.wikipedia.org/wiki/FASTQ_format)  
 [grep](http://man7.org/linux/man-pages/man1/grep.1.html)  
 [wc -l](https://ss64.com/bash/wc.html)
 </p>

 </details>

 <br>
 <br>


### Task 3:

Generate GVCF files from your samples, using the Software and steps from Day 1&2

**Alternatively:**  
 Generate GVCF files from many samples using a [bash-script with for-loops, variables and string manipulation](./bash_scripting_day1.md)

**NB:**
 - you will have to omit the BSQR step or modify it, since you do not have a validated dataset of known sites, like you would have for a model organism.
 - you need to run HaplotypeCaller with the ```-ploidy 1``` option.

  - **optional**: have a look at the [HaplotypeCaller API](https://gatk.broadinstitute.org/hc/en-us/articles/360036712151-HaplotypeCaller) ([API](https://en.wikipedia.org/wiki/Application_programming_interface) = **A**pplication **P**rogramming **I**nterface).
    Here you can find out what options are available and what they do. You dont need to know or understand any of them for now, but keep in mind that there are a lot of gears/choices hidden behind the default options.



### Task 4:

Copy the generated GVCF files into the folder folder below!  

```/proj/g2021009/private/computer_practicals/NGS_workflow_day3_4/data/GVCF```

We will run joint variant-calling on it, to generate the VCF file that we will use during the next session, or use a VCF file that we have previously generated.




### Wrap up day 1

- #### What are the main take-away messages?
- #### Which of these are relevant to the exam?
- #### How does this relate to the real world? (and our research?)

---
## Day 2

### Objectives
- Familiarise yourself with the VCF file-format
- Discuss the variant data by plotting and clustering the genotype matrix.
- Familiarise yourself with BCFtools/htslib for filtering, modifying and analysing Variant data.
- Investigate between-group allele-frequency-differences to identify candidate regions.  

the reservation code for today is:
```
uppmax2023-2-24_8
```


### Data day 2
Today we will work with the VCF file containing the samples that you processed yesterday.

### Task 1
- Copy the VCF-file from the folder below into your working directory.
```
/proj/g2020004/private/computer_practicals/NGS_workflow_day3_4/data/Mt_h37rv.vcf
```


- Today you will need to load just two libraries: bcftools and the bioinfo-tools module.

```bash

module load bioinfo-tools
module load bcftools

```

- Inspect the VCF file manually. you can use ```cat```, ```less```, ```head``` and ```grep``` for this.
  - If you're not sure what any of these do, you can read the [man page](https://en.wikipedia.org/wiki/Man_page) for these tools or google them.
  - You can also have a look at the [official specifications for vcf-fileformat 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf).
- *Can you spot a general structure?*
- *What do columns and rows represent?*
- *Can you tell how many samples there are?*


<details>
			<summary>tips</summary>
			<p>
        Dirty: (and fast to write, does not need any tool installed): grep the header-line containing the sample-names, count them.  <br>
        Clean: there's a BCFtools functionality that outputs a list of sample-names. e.g. bcftools query -l file.vcf | wc -l <br>
		What is grep doing that make it a really bad idea to use on real data?
		</p>

   </details>



### Task2
Have a look at the VCF as a graphical representation/heatmap.
Since it is bit tricky to plot and display figures on Rackham when you are just starting out, we have created the figure and added it below. If you are curious, you can look at the code that created it in detail [here](UU_NGS_heatmap.html). (It is annotated, but in Python).

![heatmap](figures/heatmap_full.png)





**Questions:**
 - _What do the X and Y axis represent?_
 - _Why did we only put "beginning" and "end" as X-axis labels?_
 - _What do the colours mean? can you guess from the context?_
 - _Can you identify duplicate samples?_
 - _How many samples do you need to identify a causative region?_
 - _Can you already spot something interesting?_


### Task3
For the third task, we are going to extract the sample names corresponding with their phenotype, using a small bashscript and BCFtools. this will later help us compare the two groups.


**Questions:**

- _Why are we comparing specific sub-phenotypes, and not just simple/complicated disease progression ?_


**tasks:**
- Extract "wildtype" and "resistant" samples-lists for a category/treatment of your choice from the table using the supplied bash script ```get_samples.sh```. Look at it here, or using ```cat``` or  ```less``` to figure out what input it needs.

its embedded below, but you can also find it at ```/proj/g2021009/private/computer_practicals/NGS_workflow_day3_4/scripts/get_samples.bash```

 ```bash  
 #! /usr/bin/env bash


 # The first argument of this command is the name of the final VCF containing all the individuals
 VCF_FILE=$1

 # The second argument of this command, $2 is the motif in the sample name that tell if it is resistant or not for the treatement you want to look at
 RES=$2

 # This command read the samples names in $VCF_FILE pass them (using |) to grep that print only the sample that contains the motif $RES
 # The names returned are then directed (with >) to be written in the file named  ${RES}.resistant.sampleList.txt (with ${RES} replaced in the real name by the second argument)
 bcftools query -l $VCF_FILE | grep $RES > ${RES}.resistant.sampleList.txt
 # The only difference with the previous one is -v is grep that do a reverse grep, taking the lines that do not have the motif $RES
 bcftools query -l $VCF_FILE | grep -v $RES > ${RES}.wild.sampleList.txt

 # now, add the population description after the sample-name
 sed -e 's/$/;resistant/' -i.bak ${RES}.resistant.sampleList.txt
 sed -e 's/$/;wild/' -i.bak ${RES}.wild.sampleList.txt

```


### Task4
For the fourth task, we will then look at the divergence in allele-frequency for each variant between groups to identify interesting variants. The metric we will be using for this is the [Fixation index](https://en.wikipedia.org/wiki/Fixation_index), though in this case, you could also see the results if you were just looking at the difference in allele-frequencies for the two groups.  
 Usually this would be done in Python, R or command line tools such as [VCFtools](https://vcftools.github.io/man_latest.html) or [plink](https://www.cog-genomics.org/plink/1.9/basic_stats). For convenience's sake, we are going to use [SNiPlay](https://sniplay.southgreen.fr/cgi-bin/analysis_v3.cgi), a webtool with a graphical user interface that also allows us to plot the results.

for this you will need the VCF, the reference and the two lists from Task3.

**Questions:**
- _Do you know what is meant with allele-frequency?_
- _Do you have an idea why we look at allele-frequency divergence and not presence/absence of variants between groups ?_
- _What can a high or low allele-frequency difference between two groups mean in this case?_
- _What could it also mean?_

**Tasks:**
- go to the SNiPlay website [here](https://sniplay.southgreen.fr/cgi-bin/analysis_v3.cgi).
- Upload the VCF file and the Reference fasta. to do this, you can download them either from the server with e.g [scp](https://linux.die.net/man/1/scp), or easier, from [here](https://drive.google.com/drive/folders/1lelvj0N75I1s853mvTz51GZDNpTQVscm?usp=sharing). If the VCF file is larger than 200 Mb, you might want to use [gzip](https://linux.die.net/man/1/gzip) to compress it before uploading.
- use the green button to upload both the ```.fa ``` reference file and the gzipped vcf (```.vcf.gz```). press the blue button to upload the files.
- go to the bottom of the page and press "submit".

<details><summary>step A</summary>

![](figures/stepa.png)

</details>


- Once you have reached the second page,  shift all individuals and the chromosome to the right, using the ```>>``` buttons.

- then, expand the ```assign individuals to populations``` section, copy and paste the content of your two lists in there. press ```submit```.


<details><summary>step B</summary>

![](figures/stepb.png)

</details>

then, click the field named ```diversity analysis```, and in the subsequent page, click ```submit```.

<details><summary>step C</summary>

![](figures/stepc.png)

</details>  

- Change the metric to 'Fst by marker', and have a look at both the figure and the table below. can you find an outlier? what does that outlier indicate?

<details><summary>step D</summary>

![](figures/stepd.png)

</details>

### Task5

**Questions:**
- Having found a potential region of interest, how would you investigate it for functional connections to the phenotype at hand, using the tools you learned about in this module or before?
- Speculate: What event did you think lead to the acquisition of resistance?

**Tasks:**
Check your hypothesis!

<details><summary>tips</summary>


You can use IGV to investigate the gene annotations for a given region. For this, you need to open IGV, then load the reference genome and the annotation file ( ```.gff```).  
Alternatively,  you can annotate the vcf using SnpEff, then look at the position with grep (section 6 and 7 of your previous practical).

</details>


### Discussion and Wrap-up.
- #### Q&A
- #### What are the main take-away messages?
- #### which of these are relevant to the exam?
