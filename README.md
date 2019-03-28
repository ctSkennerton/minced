# MinCED - Mining CRISPRs in Environmental Datasets


MinCED is a program to find Clustered Regularly Interspaced Short Palindromic
Repeats (CRISPRs) in full genomes or environmental datasets such as assembled
contigs from metagenomes. Iff you want to identify CRISPRs in raw short read 
data, in the size range of 100-200bp try using Crass (https://github.com/ctskennerton/Crass)
MinCED runs from the command-line and was derived from CRT (http://www.room220.com/crt/):
  
    Charles Bland et al., CRISPR Recognition Tool (CRT): a tool for automatic
    detection of clustered regularly interspaced palindromic repeats, BMC
    Bioinformatics 8, no. 1 (2007): 209.


## INSTALLATION

You need to install these dependencies first:
  * Java (http://www.java.com/en/download/)

there is a Makefile in the source directory so installation should be as simple as:

    cd <download_folder>
    make

To run MinCED:

    ./minced [options] file.fa

The help page can be obtained by typing:

    ./minced --help

You can get the MinCED version this way:

    ./minced --version

**NOTE: Always keep `minced` and `minced.jar` in the same folder!**

## EXAMPLES

Finding CRISPRs in the E. coli genome:

    ./minced ecoli.fna

To find repeats in short sequences, we need to decrease the minimum number of
repeats to find. For example, in 100 bp reads, we could not possibly find more
than 2 repeats:

    minced -minNR 2 metagenome.fna

The output can be large, so save it in a file:

    minced -minNR 2 metagenome.fna metagenome.crisprs

You can also save both the table output and the gff output at the same
time:

    minced ecoli.fna out.txt out.gff

## COPYRIGHT AND LICENSE

```
Copyright 2011      Florent ANGLY     <florent.angly@gmail.com>
          2013-2019 Connor SKENNERTON <c.skennerton@gmail.com>

MinCED is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
MinCED is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with MinCED.  If not, see <http://www.gnu.org/licenses/>.
```

## BUGS

All complex software has bugs lurking in it, and this program is no exception.
If you find a bug please post an issue on github https://github.com/ctSkennerton/minced/issues
