# Variant Identification Application (VIA)

VIA automates the process of variant identification, analyzing family data and returning possible candidates based on a range of specific filters and models of inheritance.

## Application Pre-Requisites
- Data Pre-Requisites
  - cleaned data file (output of Mendelian_filtering_WORK.Rmd)
  - pedigree file
- System Pre-Requisites
  - [git](https://git-scm.com/downloads)
  - [python](https://www.python.org/)

## User Guide

To install the application, clone this GitHub repository on your machine. If you have git installed on your computer, you can do this with the following command:

```
git clone https://github.com/anthonyozerov/variant-filtering.git
```

The complete application can be run through main.py with the following command:

```
python main.py
```

For greater flexibility, there are also the following optional arguments:
- **_--pedfile_ OR _-p_** : specify the absolute or relative path to the pedigree file. If no argument is specified, the application will look for a file named _Test_Ped.txt_ in the repository's directory
- **_--data_  OR _-d_** : specify the absolute or relative path to the cleaned data file. If no argument is specified, the application will look for a file named _Test_cleaned.txt_ in the repository's directory
- **_--output_ OR _-o_** : specify the file name (including extension and (optionally) the file path) for the ouput file. If no argument is specified, the application will name the output file _filtered.csv_ and place it in the repository's directory
- **_--family_ OR _-f_** : specify a certain family to output an individual csv file for. The default behaviour is to produce a single output file with variants for all families.

Any combination of these arguments can be used, and they can be chained together. For example, using all four would look like:

```
python main.py -p <file path> -d <file path> -o <file path> -f <family name>
```

## Developer Guide

### The Models (models.py)

VIA identifies variants corresponding to four models of inheritance:

#### Model 1: Homozygotes/Hemizygotes (Autosomal Recessive and X-Linked)

_TODO_

#### Model 2: Compound Heterozygotes

_TODO_

#### Model 3: De Novo

Identifies variants in affected children/siblings which are not present in either parent. Variants with an allele frequency greater than 0.0005 or a low allelic depth are removed.

#### Model 4: Autosomal Dominant

Identifies variants in affected children/siblings when _at least one_ of the parents are also affected. Variants with an allele frequency greater than 0.0005 or a low allelic depth are removed.

### Custom Classes (Family.py)

#### Person

Each sampled individual is a member of the person class. Their characteristic attributes are:
- ID (i.e. FIN5.1)
- sex
- phenotype (unaffected or affected)

#### Family

Each family belongs to the family class. Their characteristic attributes are:
- ID (i.e. FIN5)
- people - a list of Person objects 
- siblings - a list of Person objects whose relationship in the family (relative to the affected individual being studied) is 'sibling'
- mother - a Person object for the mother of the family
- father - a Person object for the father of the family
- child - a Person object for the affected individual being studied

### Custom Filters (filters.py)

_TODO_

## Change Log

### 0.1.x
- 0.1.0: First release (06.11.2021)
