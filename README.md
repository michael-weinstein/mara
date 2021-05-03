# aMino Acid Reversing App

## How to use it

The usage is pretty trivial right now, just pass your mutation as the command line argument.  At this point, only SNVs are supported, but I will be adding support for PNVs in the future.

```bash
python3 mara.py S:D614G
```

## Outputs
The output will go to STDOUT by default, but can be teed to a file or stream as needed.  An example of the output is provided in the S_D614G.json file.  The most likely codon change(s) to explain the reported amino acid will be found in the output under the key of **"best"** while less likely changes (as defined by the number of nucleotide changes required to create the observed amino acid) will be found under the key of **"other"**.

For any codon change described, the following information will be given:
* Reference and alternate codons
* Reference and alternate amino acids
* Codon start and codon end position within    ***SEE BELOW***
* Genomic positions that are mutated for the described change
* Positions within the codon that are changed as an array of booleans


### GENOMIC POSITION COORDINATE NOTES

The genomic positions given are formatted as they are given in biopython.  This means:

* Positions are given base-zero. For many variant call formats, a base-one standard is used and you will need to add 1 to the position to convert them
* The codon start and end are given using the Python standard formatting of being beginning inclusive and end non-inclusive.  A codon containing bases 9, 10, and 11 will have a start of 9 and an end of 12.
* Again, these are formatted base-zero end-noninclusive. To convert them to a base-one end-inclusive format as is used in some other applications, simply add 1 to the start position (the end will be unchanged as the change in base-value and end-inclusivity cancel one another).