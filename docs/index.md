
### Overview
### objectives
## toc

#### background & motivation
Today we will work with sequencing data derived from [_Mycobacterium tuberculosis_](https://en.wikipedia.org/wiki/Mycobacterium_tuberculosis).
 _Mycobacterium tuberculosis_ is a pathogenic bacteria, and the causative agent of [tuberculosis](https://en.wikipedia.org/wiki/Tuberculosis), a bacterial disease that kills about [1.2 million people annually](https://www.who.int/tb/publications/global_report/en/), making it the leading cause of death among infectious diseases.
 Since the treatment of tuberculosis requires sustained use of antibiotics, the rise of  multiple drug-resistant tuberculosis (MDR-TB) and extensively drug-resistant tuberculosis (XDR-TB) provides a serious concern for world health.

We will work with several "patient derived samples": strains that have been either sampled from patients where the disease responded well to antibiotic treatment, or from patients with a more problematic disease progression, where the infection did not respond to first-line antibiotics.

**Task1:**  
SSH into Rackham and request a interactive session, like you did yesterday.
The reservation code for today is ```RESERVATION_CODE ```
<details><summary>HELP</summary>
<p>

```bash
salloc -A g2019015 -t 04:00:00 -p core -n 5 --no-shell --reservation=g2019015_3 \
-M snowy &
## find your node:
squeue -u <username> -M snowy
## connect to your node:
ssh -Y <nodename>
```
</p>
</details>


<br>
<br>


We have prepared a directory with **.fastq** files from **100** samples.
each sample consists of two files, containing the forward and reverse reads.
they follow the naming scheme   ```TYPE_ID_READTYPE.FILEENDING```.
for example, sample **1337**, which has a **c**omplicated disease progression, will look like this:

```bash
C_1337_R1.fastq
C_1337_R2.fastq
```
**Task2:**
Pick one simple (S) and one complicated (C) disease progression sample and [copy]() them to your home directory.

Once you're done, please ponder the following questions. We will discuss them together in a few minutes

 - _How does this data differ from the one that you worked on before?_
 -  _Do you think you will be able to use the GATK BSQR tool on this dataset?_
 - _Can you figure out how many reads there are in one file?_
<details><summary>tips</summary>
    <p>

 [ the fastq format](https://en.wikipedia.org/wiki/FASTQ_format)  
 [grep](http://man7.org/linux/man-pages/man1/grep.1.html)  

 </p>

 </details>

 <br>
 <br>















make them look at the [HaplotypeCaller API](https://gatk.broadinstitute.org/hc/en-us/articles/360036712151-HaplotypeCaller). they dont need to read it, just realise there's a lot under the hood.
