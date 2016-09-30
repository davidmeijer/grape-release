## GARLIC  

GRAPE (Generalized Retrobiosynthetic Assembly Prediction Engine) is a software package for the retrobiosynthesis of natural small molecules, particularly those that are bacterially derived.

## Usage

PRISM is a Java package, which can be compiled then run from the command line. 

A sample command line run might look like: 

```
$ java -jar grape.jar -f smiles_file.smi -img -json
```

This runs GRAPE for a file of SMILES (a tab delimited file with names \t smiles). It will then retrobiosynthetically deconstruct these molecules and return the pieces as both an image file (png) and a json file.

To see a detailed description of all command line options, do:

```
$ java -jar grape.jar -h 
``` 

## Dependencies

The dependant libraries are located in this repository. The neo4j library will soon be replaced as currently only the json constructor (ajax) is used in the release version of GRAPE and there are no connections to a neo4j database.

GRAPE requires the CDK >= 1.5.9

## Citing GRAPE

Please cite: 

> XXXXXXXXX
