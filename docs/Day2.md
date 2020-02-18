
### objectives
- familiarise yourself with the VCF file-format
- Discuss the variant data by plotting and clustering the genotype matrix.
- familiarise yourself with BCFtools/htslib for filtering, modifying and analysing Variant data.
- investigate between-group allele-frequency-differences to identify candidate regions.
### data
Today we will work with the VCF file containing the samples that you processed yesterday.

### task1
- copy the VCF-file from the folder below into your working directory.
```
path/to/project/subproject/VCF.vcf
```
- in the same folder is a subfolder called ```scripts``` containing some premade scripts for you to use. copy the whole folder into your project directory. you can do this the same way as the vcf-file, except you will need to add the ```-r```(recursive) flag to your ```cp```-command, in order to also copy the folders contents.

- inspect the file manually. you can use ```cat```, ```less```, ```head``` and ```grep``` for this.
  - If you're not sure what any of these do, make sure to read the [man page](https://en.wikipedia.org/wiki/Man_page) for these tools.
  - You can also have a look at the [official specifications for vcf-fileformat 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf).
- *can you spot a general structure?*
- *what do columns and rows represent?*
- *can you tell how many samples there are?*
  - <details><summary>tips</summary>

      <p>
        **dirty and fast** : grep the header-line containing the sample-names, count them.  <br>
        **clean**: there's a BCFtools functionality that outputs a list of sample-names. e.g. ```bcftools query -l file.bcf | wc -l```

      </p>
   </details>



### Task2
Plot the variants as a heatmap and/or clustermap.
we have made a small script that does this for you, called ```plot_variants.py```that you can find in the folder ```/path/to/script/``` and copy over to your directory. this file is a small python script that you can run just like any bash-script. it has a few options that you can see when looking for the help-message:
you can look at the script in detail [here]()
```bash
 [in]: python plot_variants.py --help
[out]: #TODO
```
**Questions:**
 - _What do the X and Y axis represent?_
 - _What do the colours mean?_
 - _Can you identify duplicate samples?_
 - _How many samples do you need to identify a causative region?_
 - _Can you already spot something interesting?_


### Task3
For the Third task, we are going to split the VCF into multiple groups corresponding with their phenotype, using BCFtools.


**Questions:**

- _Why are we comparing specific sub-phenotypes, and not just simple/complicated disease progression ?_


**tasks:**
- Extract "wildtype" and "resistant" samples-lists for a category/treatment of your choice from the table using the supplied bash script ```get_samples.sh```. Look at it using ```cat``` or  ```less``` to figure out what input it needs.  
- split the vcf file into two files, using bcftools and the generated lists.

### Task4
for the fourth task, we will then look at the difference in allele-frequency for each variant between groups to identify interesting variants.

**Questions:**
- _Do you know what is meant with allele-frequency?_
- _Do you have an idea why we look at allele-frequency difference and not presence/absence of variants between groups ?_
- _What does a high or low allele-frequency difference between two groups mean in this case?_
- _What could it also mean?_

**Tasks:**
for both files, extract allele-frequencies using BCFtools.
compute and plot the allele-frequency delta using the provided script ```delta_af.py```

**Questions:**



### Discussion and Wrap-up.
