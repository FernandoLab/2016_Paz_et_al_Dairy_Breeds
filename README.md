To recreate the anlaysis used in the Paz et al. manuscript, “Rumen microbial composition in Holstein and Jersey cows is different under same dietary condition and is not affected by sampling method”, there are two steps (follow the guidelines below). All of the commands to generate the manuscript outputs have been ran on Mac OS X 10.10.5 (others systems should work) with 8 GB RAM. No root access is needed. Analyses should also work in a linux environment as well if the linux versions of USEARCH and [anaconda package manager](https://www.continuum.io/downloads) are used. The two dependencies needed are X11 (remember if logging onto a server) and perl (version shouldn't matter).
 
  1. Run the bash script to create a virtual enironment and download/install programs **LOCALLY** with the anaconda package manager. This will recreate the same enivronment used during the analyses of the data in the manuscript.
  
  2. Render the R Markdown file with knitR to recreate the workflow and outputs.

Due to licensing constraints, USEARCH could not be included in the setup. To obtain a download link, go to the USEARCH [download page](http://www.drive5.com/usearch/download.html) and select version USEARCH v8.1.1756 for Mac OSX. **A link (expires after 30 days) will be sent to the provided email. Use the link as an argument for shell script below.**

Simply download the bash script from the github repository and run it (provide the link to download your licensed USEARCH version as an argument for setup.sh):

  1. wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/setup.sh
  2. chmod 775 setup.sh
  3. ./setup.sh usearch_link

**Anaconda is downloaded and prompts you during installataion of the packages above. The prompts are as follows:**

  1. Press enter to view the license agreement
  2. Press enter to read the license and q to exit
  3. Accept the terms
  4. Prompts you where to install anaconda. Simply type anaconda to create a directory within the current directory. Should be: [/Users/user/anaconda2] >>> anaconda
  5. No to prepend anaconda to your path. Choosing yes does not impact the installation.
  6. Will be asked a few times if you wish to proceed with installing the packages, agree to it.
  7. After installation, enter ‘source anaconda/bin/activate projectEnv’ on the command line to activate the virtual enviornment with all dependencies.

To convert the R markdown to html use the command: **render(“dairy_breeds.Rmd”)**. To start a R session and run the workflow, use these commands from within the direcotry you initiated installation:

 1. source anaconda/bin/activate projectEnv
 2. R
 3. install.packages("rmarkdown", repos='http://cran.us.r-project.org')
 4. install.packages("knitr", repos='http://cran.us.r-project.org')
 5. source("http://bioconductor.org/biocLite.R")
 6. biocLite("Heatplus", ask=FALSE, suppressUpdates=TRUE)
 7. library(rmarkdown)
 8. library(knitr)
 9. render("dairy_breeds/dairy_breeds.Rmd")

The rendered html version can be found [here](./dairy_breeds.html)
