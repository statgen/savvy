# M3VCF Specification

## Variable Length String (VLS) Encoding
```
+----------+~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
| SSSSSSSS | AAAAAAAA |         STRING_DATA         |
+----------+~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
* S: Size of string byte.
* A: Additional size of string bytes. All size bytes sum to total string size. An additional size byte can only be preceded by a size byte of 255.
* STRING_DATA: String payload stored in S bytes.
```


## Header Format
```
+-----------------------------------------------------------------------------------+
|           "m3vcf" string in binary + 4 bytes for version (Major.minor)            |
+-----------------------------------------------------------------------------------+
| 01101101 00110011 01110110 01100011 01100110  MMMMMMMM MMMMMMMM mmmmmmmm mmmmmmmm |
+-----------------------------------------------------------------------------------+


+vvvvvvvvv+-------------------------------------+-------------------------------------+VVVVVVVVVVVVVVVVVVVVVVVVVVVVV+
|  CHROM  | PPPPPPPP PPPPPPPP PPPPPPPP PPPPPPPP | SSSSSSSS SSSSSSSS SSSSSSSS SSSSSSSS |     SAMPLE_ID_ARRAY ...     |
+vvvvvvvvv+-------------------------------------+-------------------------------------+VVVVVVVVVVVVVVVVVVVVVVVVVVVVV+
* CHROM: Chromosome string stored has VLS.
* P: Ploidy level stored in 1 byte.
* S: Number of samples stored in 4 bytes.
* SAMPLE_ID_ARRAY: Array of length S that stores sample ID's in VLS encoding.

+----------+VVVVVVVVVVVVVVVVVVVVVVVVVVVVV+
| SSSSSSSS |     META_FIELDS_ARRAY ...   |
+----------+VVVVVVVVVVVVVVVVVVVVVVVVVVVVV+
* S: Number of metadata fields.
* META_FIELDS_ARRAY: Array of length META_FIELDS_CNT that stores metadata fiels in VLS encoding.
```


## Block Format
```
+-------------------------------------+---------------------------------------+VVVVVVVVVVVVVVVVVVVV+
| MMMMMMMM MMMMMMMM MMMMMMMM MMMMMMMM | NNNNNNNN  NNNNNNNN  NNNNNNNN NNNNNNNN | SAMPLE_MAPPING ... |
+-------------------------------------+---------------------------------------+VVVVVVVVVVVVVVVVVVVV+
+-------------------------------------------------------------------------+vvvvvvvvvv+vvvvvvvvv+vvvvvvvvv+VVVVVVVVVVVVVVVVVVVVVV+VVVVVVVVVVVVVVVVVVVV+
| PPPPPPPP PPPPPPPP PPPPPPPP PPPPPPPP PPPPPPPP PPPPPPPP PPPPPPPP PPPPPPPP |  MARKID  |   REF   |   ALT   | META_VALUE_ARRAY ... | UNIQUE_HAP_ROW ... |
+-------------------------------------------------------------------------+vvvvvvvvvv+vvvvvvvvv+vvvvvvvvv+VVVVVVVVVVVVVVVVVVVVVV+VVVVVVVVVVVVVVVVVVVV+
...


* M: Number of rows in unique haplotype matrix (number of markers in block).
* N: Number of columns in unique haplotype matrix.
* SAMPLE_MAPPING: Array of integers with a byte width of ceil(log2(N + 1) / 8). Array contains sample count times ploidy level elements.
* P: Chromosome pos stored in 8 bytes.
* MARKID: Marker ID string stored has VLS.
* REF: Reference haplotype stored has VLS.
* ALT: Alternate haplotype stored has VLS.
* META_VALUE_ARRAY: Array of VLS encoded strings that stores metadata for each marker. Values correspond to META_FIELDS_ARRAY in header.
* UNIQUE_HAP_ROW: Array of '.', '0' or '1' bytes of size N.

```