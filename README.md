# Generating a subset of REAL that is related to the DSIP library

The query.py is a simple Python script that does a bulk similarity search of
one set of molecules against another. 

The DSIP molecules are read, and their fingerprints stored in an array.
The molecules to filter are then read one-by-one and RDkit's
[BulkTverskySimilarity](http://rdkit.org/docs/source/rdkit.DataStructs.cDataStructs.html?highlight=bulktverskysimilarity#rdkit.DataStructs.cDataStructs.BulkTverskySimilarity)
function is used to generate the similarity scores for that molecule against all the DSIP 
molecules and then retain only those shose best score is within the target range. 

We use an asymmetric function so that we can find larger molecules that match just part
of the target molcule.

Run the script with something like this:
```
python query.py dsip+direct_frags.smi targets.sdf 0.795 0.805 0.0
```
The parameters:
1. the filename of the DSIP molecules (.smi or .sdf)
2. the filename of the queries (REAL molecule as .smi or .sdf)
3. the lower cuttoff for the similarity score
4. the upper cuttoff for the similarity score
5. the value for alpha (beta is calculated as 1.0 - alpha)

In order to define an optimal set the following parameters need consideration:

1. The DSIP molecules to consider. Just the molecules themselves or some of their fragments.
2. The fingerprint to use.
3. The alpha and beta values for the Tversky algorithm
4. The similarity threshold.


## Investigating the REAL building blocks

The first study is to try to filter the REAL building blocks and find ones that are related to the DSIP library.

The first attempt used the DSIP molecules and all child fragments (recursively) that either had a ring or
had 5 or more heavy atoms. This resulted in every buildign block matching. Clearly we generate so may promiscuous
fragments that at least one of them is present in every molecule. We need to be more specific.

Two approaches were then used:

1. Using the DSIP molecules and no fragments
2. Using the DSIP molecules and direct child fragments that had either:
2.1. a ring and at least 2 non-ring atoms
2.2. no ring but at least 5 heavy atoms.

To assess the parameters searches were run with  values for alpha (beta was calculated as 1.0 -alpha) of
0.0 (completely asymmetric) and 0.2 (partly asymmetric). For each of these slices of the results around
a given similarity score (+- 0.05) so that the types of molecules at different similary scores could be visualised.

The results are in [slices](). Those with names starting with `dsip_` are with just the DSIP library, those 
starting with `dsip+direct_frags_` are for the DSIP plus direct child fragments as described above.
The next part of the filename is `a0.0` or `a0.2` which defines the value fo rthe alpha parameter.
The final part is like `t0.7` where the 0.7 bit is the simlarity score used (e.g. in that case
molecules with scores between 0.695 and 0.705.

The output first lists the parameters. e.g something like:
```
Target: dsip_mols.smi Query: 2020q3-4_Enamine_REAL_reagents_SDF.sdf MinThreshold: 0.595 MaxThreshold: 0.605 Alpha: 0.0
Fingerprinting targets took: 0.0659940242767334
Number of targets: 768
```

Then has pairs of lines like this:
```
Nc1nc(-c2ccccc2)cs1 0.6
CC(C)Oc1ccc(-c2csc(N)n2)cc1
```

Where the first line is the molecule being compared (e.g. from the REAL buidling blocks) followed by the 
best similarity score foudn for that molecule, and the second line is that molecule or fragment from DSIP
that generated that best score. It is formatted in this way so that it can be pased into
[CDK Depict](https://www.simolecule.com/cdkdepict/depict.html) so that you can see the molecule being 
considered and it's most similar molecule in DSIP.

At the end of the file is some summary statistics e.g.
```
Queries took: 32.88010025024414
Found 2108 molecules
```




