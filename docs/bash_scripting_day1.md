


### Basic for-loop:


```bash
for i in iterable; #e.g. file in folder: for i in ./*.fileending
do                     # beginning of loop
  echo "do something"; # the action that should be performed for each i in iterable
                       # in the for-loop, i can be accessed using the ${i} variable.
done                   # end of loop
```

### command-line variables in bash:

numerical bash variables correspond to positional command line arguments:

consider the following script,
```cheese.sh```:

```bash
#!/bin/bash

echo ${1} # print the first positional command

if [[ ${2}="cheese" ]]; then  # if the second positional command is "cheese"
  echo "cheddar"              # print "cheddar"
fi

```

running this script as-is will not result in anything being printed to the output.
but adding another word after the script, will result in this word being printed, the first position after the script corresponding to
```$1```.  ( or ```${1}```):


```sh
 [in]: bash ./cheese.sh mouse
[out]: mouse
```

running :

```sh
 [in]: bash ./cheese.sh mouse apple
[out]: mouse
```

will not result in "apple" being printed, since the second positional argument (```$2```) is wrapped into a basic if-clause: only if the second word is "cheese", the script will print the word "cheddar":

```sh
 [in]: bash ./cheese.sh mouse cheese
[out]: mouse
[out]: cheddar
```

### String manipulation:

```bash
 [in]: test_variable="path/to/SampleA_R1.fastq" # define variable
 [in]: echo ${test_variable}                    # print variable
[out]: path/to/SampleA_R1.fastq
```        
the most important string manipulation is probably **substring removal**:<br>
```%substring``` removes the **shortest matching substring** from the **back** of the string.

```bash
 [in]: echo ${test_variable%.fastq}  #remove the fileending ".fastq"           
[out]: path/to/SampleA_R1
 [in]: echo ${test_variable%.fastq}.sam  #you can add new strings if you do so outside of the curly brackets           
[out]: path/to/SampleA_R1.sam
                  
```

```#substring``` removes the **shortest matching substring** from the **front** of the string.

 so in our case, ```*/```  could match either ```path/``` or ```path/to/```, but ```#substring``` only matches the shortest possible version. ( the star is a [regular expression](https://en.wikipedia.org/wiki/Regular_expression))


```##substring``` and ```%%substring``` match and remove the **longest** possible substring.


```bash
 [in]: echo ${test_variable##*/}  #remove the filepath           
[out]: SampleA_R1.fastq  
 [in]: echo totally/new/path/to/${test_variable##*/}  #you can add new strings if you do so outside of the curly brackets           
[out]: totally/new/path/to/SampleA_R1.fastq

```
