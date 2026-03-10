# Finestructure Version 4

## About fineSTRUCTURE

fineSTRUCTURE is a fast and powerful algorithm for identifying population structure using dense sequencing data.  By using the output of ChromoPainter as a (nearly) sufficient summary statistic, it is able to perform model-based Bayesian clustering on large datasets, including full resequencing data, and can handle up to 1000s of individuals. Full assignment uncertainty is given.

finestructure works on Linux and Mac. You can also compile it for Windows if you set up the required toolchains, but you need to know what you are doing and be comfortable with the command line.

A Stochastic optimization routine is available for performing faster EDA and dealing with larger datasets - see FAQ under "What if my dataset is too big for MCMC".

Important Note: fs4 includes ChromoPainter, which have different licences and authors. Both are free for Academic use only, and explicitly exclude commercial applications. See the file COPYING in the download for details.

Please consider [Registering](https://forms.office.com/Pages/ResponsePage.aspx?id=MH_ksn3NTkql2rGM8aQVGwBpDrkt7zVLlSvbqowMvq1UQzNXV1hMQkE1QzBNTDJPS0JEWFE4TDFaVi4u), if you have not already done so.

This codebase can compile to produce `fs`, the binary executable for finestructure (and chromopainter), as well as `mixPainter` which is a limited version made for simple use of [badMIXTURE](https://github.com/danjlawson/badMIXTURE).

## Technical details:

NB Installation instructions are below: Linux and Mac OS X should work fully.

### LICENCE

* The "fs" code is for non-commercial purposes only.
* It is free to use for Academic, personal and non-profit purposes.
* Attribution to software creator must be made:
    	* Acknowledgement in personal or other non-commercial work. 
		* For academic use, citation of the appropriate article(s), currently: [Lawson, Hellenthal, Myers & Falush 2012, PLoS Genetics e1002453 "Inference of population structure using dense haplotype data"](https://journals.plos.org/plosgenetics/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1002453).
* For commercial licencing, contact the author dan.lawson <at> bristol.ac.uk.
* We are likely to agree to modification and distribution for non-commercial use, but please ask.
* fs makes use of the non-free software ChromoPainter version 1, which has been modified with permission of the author Garrett Hellenthal (ghellenthal <at> gmail.com). You are bound by the terms of that software SEPARATELY. It uses the same licence as fs but the copyright belongs to Garrett Hellenthal and permission to use the software outside of these terms must be arranged with both authors separately.

See LICENCE for further details.

### REQUIREMENTS

The requirements are a "recent" version of the GCC compiler, for the GLIBC libraries. These are typically present by default.

There are probably other requirements which are present on most systems by default.  If you encounter any problems let us know. 

### INSTALLATION

A binary for linux and mac is provided, which you may just be able to use. However some users have different versions of some core C++ libraries that make compilation desirable.

fs4 comes with an **installation script** `fs_install.sh`. This is optional but will help novices get the correct version for their computer as well as adding it to their PATH variable.

### COMPILATION 

If you have a similar enough build-chain to me you can use:

```{sh}
./configure
make
make install #optional
```

If you need to reconfigure the configure options, to change the 

If you are compiling on a Mac, the default C compiler (clang) does not support multi-threading.  If you would like to use another C compiler that is not the default, you can use:

```{sh}
./configure CXX=g++-11
```

(substitute your CXX version). This is implemented in ./configure.mac so that you don't have to remember each time.

If your toolchain is a little different, try recreating the configure file:

```{sh}
./autogen.sh
```

This script requires a more complex toolchain, noteable `automake` and `autoreconf`, which can usually be installed via `brew` (for mac), apt-get (for Ubuntu), or other package manager.

To compile `mixPainter`, ince you are able to compile fs with the above `make` command, simply run:

```{sh}
make mixPainter
```

### Dependencies

Dependencies for the command line version are (from version 4.1.0) nearly non-existent!

* Linux:
  1. GCC c++ compiler (package `build-essential` in ubuntu)
  2. Automake (optional)
* Mac OS: Either:
  1. Xcode c++ "clang" compiler, installed via "Xcode Command Line Tools" with `xcode-select –install` (See [e.g. This Howto](https://mac.install.guide/commandlinetools/4.html); but note that this does not enable parallel computation outside of HPC mode.)
  2. or GCC c++ compiler (package `gcc` in [brew](https://brew.sh/)) for full functionality; to install use`brew install gcc`, or if you are starting from scratch:
```{sh}
brew update
brew upgrade
brew info gcc # Gets information about it
brew install gcc
```
	3. Automake (optional; install all required tools with `brew install autoconf automake libtool`)
 
### FURTHER INFORMATION

You need to prepare your data in CHROMOPAINTERv2 format. This is not trivial from some file formats.

Importing from:
* **Impute2** format: use the provided `impute2chromopainter.pl` script.
* **Beagle** format: use the provided `beagle2chromopainter.pl` script. Note that this is for early versions of Beagle; they now use VCF.
* **VCF** format:

This is nearly the same as PHASE format, but HAS BEEN UPDATED since chromopainter v 0.0.5. We provide some tools for this in the scripts directory.

Run "fs" for help. Examples are included in the "examples" directory; it is recommended to work through the examples to establish how to run this program on your own dataset.

IMPORTANT: If you have a small dataset, you can run 
```{sh}
fs project.cp -phasefile <data> -recombfile <recombination map file> -go
```
and it will do everything with default settings that should work. BUT IF YOU HAVE A LARGE DATA SET THIS IS GOING TO TAKE A LONG LONG TIME. You will want to parallelise the work, which is done with the "-hpc" flag. You can then run commands on your own HPC or in parallel on a multi-core machine.

The examples show you how to do this!

[www.paintmychromosomes.com](www.paintmychromosomes.com) is the place to start for help.  Visit the FAQ page for standard issues.

fs4 (finestructure) is written by Daniel Lawson (dan.lawson@bristol.ac.uk) COPYRIGHT University of Bristol 2022.

## Changelog

### 2026-03-10 — Ne/mu aggregation changed from mean to median

In `fsproject.cpp` (`combineCpEm()`), the aggregation of per-individual Ne and mutation rate (mu) estimates from ChromoPainter EM runs (stages 1 and 6) was changed from arithmetic mean to median. This makes the global parameter estimates more robust to outlier individuals whose EM runs converge to extreme values.

### 2026-03-09 — New and updated conversion scripts in `scripts/`

Three new Perl scripts were added:

* **`scripts/convertrecfile_v2.pl`** — Updated version of `convertrecfile.pl`. Fixes the hapmap-format column indices to use columns 1 and 3 (0-indexed) instead of 2 and 4, matching the 3-column linkage map format produced by SHAPEIT5.

* **`scripts/impute2chromopainter_v2.pl`** — Updated version of `impute2chromopainter.pl`. Input is now read from STDIN (pipe or redirect) rather than as a positional filename argument, enabling `zcat file.haps.gz | perl impute2chromopainter_v2.pl prefix` usage. Adds `-f` flag to prepend the fineSTRUCTURE-required `0` header line. It also uses a more efficient data structure and requires less ram.

* **`scripts/impute2chromopainter_v3.pl`** — Extends v2 with `-hap <file>` (read `.haps`/`.haps.gz` directly without piping) and `-legend <file>` (supply a bcftools `--haplegendsample` legend file for `.hap` inputs that lack SNP metadata columns). This would output .phase file with the first 5 info columns.

### 2026-03-09 — Invariant SNP filtering based on ID file

Added `filterInvariantSNPs()` (`cp/ChromoPainterData.c`) which removes SNPs that are invariant among the individuals included via `-idfile` before painting begins. The filter runs after `assignRecMap` so that recombination rates across removed sites are correctly merged as a Morgan-distance-weighted average into a single effective rate for the collapsed interval. Missing/gap alleles (codes 8 and 9) are ignored during the invariant test. A fatal error is raised if fewer than 2 SNPs remain after filtering. The call site is in `chromopainter()` (`cp/ChromoPainterMutEM.c`), between `assignRecMap` and `makeHeaders`, so output file headers reflect the filtered SNP count.

### 2026-03-09 — ChromoPainter segfault fix on truncated recombination file

Fixed a crash (segmentation fault) that occurred when ChromoPainter encountered an error during processing (e.g. a recombination file shorter than expected). The error-handling path uses `setjmp`/`longjmp`; two bugs combined to cause the crash:

1. **`DestroyIds` lacked a NULL guard** (`cp/ChromoPainterData.c`). Unlike the sibling functions `DestroyData`, `clearDonors`, and `clearCopyvec` which all guard against a NULL argument, `DestroyIds` dereferenced its pointer unconditionally. If `longjmp` fired before `Ids` had been assigned, cleanup passed NULL and the function crashed.

2. **`Ids`, `Data`, `Donors`, `Copyvec` were not `volatile`** in `chromopainter()` (`cp/ChromoPainterMutEM.c`). The C standard (§7.13.2.1) states that local variables modified between `setjmp` and `longjmp` have indeterminate values after `longjmp` unless declared `volatile`. The compiler may hold these pointer values in registers; when `longjmp` restores the register snapshot taken at `setjmp` time, the pointers revert to NULL (or stale values), so the cleanup function received bad pointers. Declaring them as `T * volatile` forces the compiler to keep them in memory, ensuring cleanup always sees the current value.

###############################

How to install perl dependencies?

#### Local install - not root privileges:
```{sh}
## NB: This requires a requirement to run perl scripts with `perl -Mlocal::lib /path/to/script.pl`
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Switch)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(VCF)'
```

#### Global install - uses root privileges:
```{sh}
## NB: This allows running perl scripts with `/path/to/script.pl`
sudo perl -MCPAN -e 'CPAN::install(Switch)'
sudo perl -MCPAN -e 'CPAN::install(VCF)'
```

```{sh}
## Get some data in VCF format
git clone git@github.com:danjlawson/pcapred.ref.git
cp pcapred.ref/inst/extdata/1000G_tinysubset.* .
gunzip 1000G_tinysubset.bim.gz
plink1.9 --bfile 1000G_tinysubset --recode vcf --out 1000G_tinysubset_unphased
## Process each chromosome separately:
for chr in `seq 1 22`; do
	## First phase the data:
	java -jar $HOME/bin/beagle.28Jun21.220.jar gt=1000G_tinysubset_unphased.vcf out=1000G_tinysubset_chr$chr chrom=$chr
	## Convert it to chromopainter format via the safe VCF route:
	gunzip 1000G_tinysubset_chr$chr.vcf.gz
	perl -Mlocal::lib ~/bin/vcf2cp.pl 1000G_tinysubset_chr$chr.vcf 1000G_tinysubset_chr$chr
	## Make a suitable recombination map:
	makeuniformrecfile.pl 1000G_tinysubset_chr$chr.phase 1000G_tinysubset_chr$chr.rec
done
## Run a combined finestructure analysis:
## NB The format {1..22} is bash specific and you may have to list the files individually.
fs 1000G_tinysubset_test.cp -phasefiles 1000G_tinysubset_chr{1..22}.phase -idfile 1000G_tinysubset_chr1.ids -recombfiles 1000G_tinysubset_chr{1..22}.rec -go
```


