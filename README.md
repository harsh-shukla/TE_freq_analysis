# TE population frequency analysis

This repository contains bash and python scripts used for calling presence/absence of TEs in a sequenced genome(illumina)

## Description of sub-folders

```
00-Processing_Data
```
This folder contain scripts used for pre-processing the data (downloading,trimming)

```
01-Mapping_and_SVCalling
```
This folder contains a bash script that contains commands used to map reads and call( &filter ) SVs.

```
02-Updating_boundaries
```
This folder contains scripts used to update TE boundaries

```
03-TE_Calling
```
This folder contains the scripts finally used to call absence/presence of TEs using the updated TE annotation(From _02-Updating\_boundaries_) and mapped bam(From _01-Mapping\_and\_SVCalling_)

## Authors

* **Harsh Shukla**

## Reference

Yuheng Huang, Harsh Shukla, Yuh Chwen G. Lee. "Species-specific chromatin landscape determines how transposable elements shape genome evolution." elife 11 (2022): e81567. [HTML](https://doi.org/10.7554/eLife.81567)

## Contact

For any queries I can be reached at : *harsh.g.shukla@gmail.com*

