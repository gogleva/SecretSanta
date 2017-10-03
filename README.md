[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](http://www.repostatus.org/badges/latest/concept.svg)](http://www.repostatus.org/#concept)
[![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Last-changedate](https://img.shields.io/badge/last%20change-2017--10--03-yellowgreen.svg)](/commits/master)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/)
[![Travis-CI Build Status](https://travis-ci.org/gogleva/SecretSanta.svg?branch=master)](ht#tps://travis-ci.org/gogleva/SecretSanta)

## 1. Background

The **SecretSanta** package provides an R interface aiding integrative prediction of extracellular proteins secreted via classical secretion pathways.

Secretome prediction often involves multiple steps. Typically, it starts with prediction of short signal peptides at the N-terminal end of a protein. Next, it is crucial to ensure the absence of motifs and domains preventing the protein from being secreted despite the presence of the signal peptide. These sequences include transmembrane domains, short ER lumen retention signals,and mitochondria/plastid targeting signals.

A number of excellent command line tools and web-interfaces exist to perform predictions of individual motifs and domains ([signalp](http://www.cbs.dtu.dk/services/SignalP/), [targetp](http://www.cbs.dtu.dk/services/TargetP/), [TMHMM](http://www.cbs.dtu.dk/services/TMHMM/), [WolfPsort](https://github.com/fmaguire/WoLFPSort)), however the interface allowing to combine the outputs in a single flexible workflow is lacking.

**SecretSanta** package attempts to bridge this gap. It provides wrapper functions around existing command line tools for prediction of signal peptides and protein subcellular localisation. The wrappers are designed to work together by producing standardized output. This allows to pipe results between individual predictors easily to create flexible custom pipelines and also to compare predictions between similar methods.

To speed-up processing of large input fasta files initial steps of the pipeline are automatically run as an embarrassingly parallel process when the number of input sequences exceeds a certain limit.

Taken together **SecretSanta** provides a platform to build automated multi-step secretome prediction pipelines that can be applied to large protein sets to facilitate comparison of secretomes across multiple species or under various conditions.

## 2. External dependencies

SecretSanta relies on a set of existing command line tools to predict secreted proteins. Please install them and configure according to the listed instructions. Due to limitations imposed by the external dependencies, some of SecretSanta wrapper functions won't work in Windows, however should be fully functional on Linux and Mac.

### 2.1 Automatic installation of external dependencies

Download the external dependencies:

-   **siganlp-2.0** <http://www.cbs.dtu.dk/cgi-bin/sw_request?signalp+2.0>
-   **signalp-3.0** <http://www.cbs.dtu.dk/cgi-bin/sw_request?signalp+3.0>
-   **signalp-4.1** <http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp>
-   **targetp-1.1** <http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?targetp>
-   **WoLFpsort** <https://github.com/fmaguire/WoLFPSort.git>
-   **tmhmm-2.0** (<http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm>)

Place all the tarballs in a dedicated directory and run the following installation script inside it: <https://gist.github.com/gogleva/9ee60038d34327e300edf082a786fa55>

### 2.2 Manual installation of external dependencies

#### Tools for prediction of signal peptides and cleavage sites:

-   **signalp-2.0**
    -   This version can run under IRIX, IRIX64, Linux, OSF1, SunOS.
    -   Download stand alone signalp2.0 <http://www.cbs.dtu.dk/cgi-bin/sw_request?signalp+2.0>
    -   Unpack the archive

    ``` sh
    tar -zxvf signalp-2.0.Linux.tar.Z
    cd signalp-2.0
    ```

    -   Edit "General settings" at the top of the **signalp** file. Set value of 'SIGNALP' variable to be path to your **signalp-2.0** directory. Other variables usually do not require changes. We will not use plotting functions from signalp, so **gnuplot**, **ppmtogif** and **ghostview** are not required. For more details please check `signalp-2.0.readme`.
    -   Since, we want to be able to run different versions of **signalp**, including the legacy ones, it is important to be able to discriminate between them. R is oblivious to shell aliases, so we will simply rename the **siganlp** script:

    ``` sh
    mv signalp siganlp2
    ```

-   **signalp-3.0**
    -   This version will run on the most common UNIX platforms.
    -   Download stand alone signalp3.0 <http://www.cbs.dtu.dk/cgi-bin/sw_request?signalp+3.0>
    -   Unpack the archive

    ``` sh
    tar -zxvf signalp-3.0.Linux.tar.Z
    cd signalp-3.0
    ```

    -   Similar to **signalp-2.0**, edit "General settings" at the top of the signalp file. Set value of 'SIGNALP' variable to be path to your **signalp-3.0** directory. Other variables usually do not require changes. For more details please check `signalp-3.0.readme`.
    -   Rename **signalp** script to avoid further confusion between the versions:

    ``` sh
    mv signalp signalp3
    ```

-   **signalp-4.1** - the most recent version
    -   This version can run under Windows, OS X (Macintosh) and Linux.
    -   Download stand alone signalp4.0 <http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp>
    -   Unpack the archive

    ``` sh
    tar -zxvf signalp-4.1.Linux.tar.Z
    cd signalp-4.1
    ```

    -   Edit "General settings" at the top of the **signalp** file. Set values for 'SIGNALP' and 'outputDir' variables. For more details please check `signalp-4.1.readme`.
    -   Rename **signalp** script to avoid further confusion between the versions:

    ``` sh
    mv signalp signalp4
    ```

    #### Tools for prediction of protein subcellular localization:

-   **taretp-1.1**
    -   **tatgetp-1.1** will run on the most common UNIX platforms
    -   Download stand alone targetp <http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?targetp>
    -   Unpack the archive:

    ``` sh
    tar -zxvf targetp-1.1b.Linux.tar.Z
    cd targetp-1.1
    ```

    -   Edit the paragraph labeled "GENERAL SETTINGS, customize" at the top of the **targetp** file. Set values for 'TARGETP' and 'TMP' variables. Ensure, that the path to **targetp** does not exceed 60 characters, otherwise **targetp-1.1** might fail.
-   **WoLFPsort**
    -   Clone WoLFPsort

    ``` sh
    git clone https://github.com/fmaguire/WoLFPSort.git
    cd WoLFPSort
    ```

    -   Copy the binaries from the appropriate platform specific binary directory `./bin/binByPlatform/binary-?` to \`./bin/\`\`
    -   For more details please check the `INSTALL` file.
    -   The most important script we need **runWolfPsortSummary** has a bulky name, we will rename it to simply **wolfpsort** for the future convenience:

    ``` sh
    mv runWolfPsortSummary wolfpsort
    ```

#### Tools for prediction of transmembrane domains

-   **tmhmm-2.0**
    -   **tmhmm-2.0** will run on the most common UNIX platforms
    -   Download stand alone tmhmm (<http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm>)
    -   Unpack the archive:

    ``` sh
    tar -zxvf tmhmm-2.0c.Linux.tar.gz
    cd tmhmm-2.0c
    ```

    -   Set correct path for Perl 5.x in the first line of `bin/tmhmm` and `bin/tmhmmformat.pl` scripts.
    -   For more details please check the `README` file.

#### Organise access to the external dependencies

The best option would be to make all the external dependencies are accessible from any location. This requires modification of `$PATH` environment variable.

To make the change permanent, edit `.profile`: Open ./profile:

``` sh
gedit ~/.profile
```

Add a line with all the path exports. In this example all the dependencies are installed in the `my_tool` directory:

``` sh
export PATH=
"/home/my_tools/signalp-4.1:\
/home/my_tools/signalp-2.0:\
/home/my_tools/signalp-3.0:\
/home/my_tools/targetp-1.1:\
/home/tmhmm-2.0c/bin:\
/home/my_tools/WoLFPSort/bin:\
$PATH"
```

Reload `/.profile`:

``` sh
. ~/.profile
```

Reboot, to make changes visible to R.

### 3. Installation

To install **SecretSanta** package:

``` r
library("devtools")
install_github("gogleva/SecretSanta", build_vignettes = TRUE)
library("SecretSanta")
browseVignettes(package = "SecretSanta") # manual
```

### Reporting bugs

email <anna.gogleva@slcu.cam.ac.uk> about bugs and strange things.
