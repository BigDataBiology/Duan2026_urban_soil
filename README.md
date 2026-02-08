# Long-read metagenomic sequencing reveals novel lineages and functional diversity in urban soil microbiome

This repository contains files and scripts to generate analysis and figures in the manuscript _Long-read metagenomic sequencing reveals novel lineages and functional diversity in urban soil microbiome_:

> Yiqian Duan, Anna Cuscó, Chengkai Zhu, Yaozhong Zhang, Alexandre Areias Castro, Xinrun Yang, Jiabao Yu, Gaofei Jiang, Xing-Ming Zhao, Luis Pedro Coelho

See also the **MAG collection** at the [Urban soil MAG collection
website](https://urban-soil-mags.netlify.app/).

## Structure

The folder `resource_generation` contains scripts to generate the MAG and gene catalog from raw data. The scripts are provided for transparency and reproducibility, but we recommend using the precomputed data (deposited at Zenodo or included here, see below).

The folder `analysis` contains pre-computed files and scripts to run the analysis and generate figures included in the manuscript. The scripts are written in Python (depenencies listed below). Generally speaking, these do not require large computational resources and interested users can run them on their own machines and adapt them to perform follow-up analyses.

## Dependencies

The following are required for the scripts (other versions may work, we list the ones that were used).

| **Software** | **Availability** |
| :---: | :---: |
| Chopper (v.0.3.0) | https://github.com/wdecoster/chopper |
| Porechop (v 0.5.0) | https://github.com/bonsai-team/Porechop_ABI |
| NGLess (v.1.5) | https://github.com/ngless-toolkit/ngless |
| Flye (v.2.9.2) | https://github.com/mikolmogorov/Flye |
| Medaka (v.1.9.1) | https://github.com/nanoporetech/medaka |
| Polypolish (v.0.5.0) | https://github.com/rrwick/Polypolish |
| MaSuRCA (v.4.1.0) | https://github.com/alekseyzimin/masurca |
| SemiBin2 (v.1.5.1) | https://github.com/BigDataBiology/SemiBin |
| Checkm2 (v.1.0.1) | https://github.com/chklovski/CheckM2 |
| GUNC (v.1.0.6) | https://github.com/grp-bork/gunc |
| Barrnap (v.0.9) | https://github.com/tseemann/barrnap |
| tRNAscan (v.2.0.12) | https://github.com/UCSC-LoweLab/tRNAscan-SE |
| CoverM (v.0.7.0) | https://github.com/wwood/CoverM |
| dRep (v.3.5.0) | https://github.com/MrOlm/drep |
| fetchMGs (v.2.0.1) | https://github.com/motu-tool/FetchMGs |
| GTDB-tk (v.2.4.1) | https://github.com/Ecogenomics/GTDBTk |
| antiSMASH (v.7.0.0) | https://github.com/antismash/antismash |
| BiG-SCAPE (v.2.0) | https://github.com/medema-group/BiG-SCAPE |
| Prodigal (v.2.6.3) | https://github.com/hyattpd/Prodigal |
| Eggnog-mapper (v.2.1.12) | https://github.com/eggnogdb/eggnog-mapper |
| RGI (v.6.0.3) | https://github.com/arpcard/rgi |
| GMSC-mapper (v.0.1.0) | https://github.com/BigDataBiology/GMSC-mapper |
| CD-HIT (v.4.8.1) | https://github.com/weizhongli/cdhit |
| geNomad (v.1.8.1) | https://github.com/apcamargo/genomad |

## Data Availability

### Database

These databases are used in the construction and analysis of the catalogue.

| **Database** | **Availability** |
| :---: | :---: |
| SPIRE | https://spire.embl.de |
| AntiFam (v.7.0) | ftp://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/ |
| The Conserved Domain Database | https://ftp.ncbi.nih.gov/pub/mmdb/cdd |
| GTDB R226 | https://gtdb.ecogenomic.org/ |

### Preprocessed data

MAG catalogue & annotations: The MAG catalogue and its annotations are available at zenodo

Preprocessed data: For convenience, the preprocessed files are available under the `anlysis/pre-calculated_data` folder.
