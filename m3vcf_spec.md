# M3VCF Specification

## Variable Length String (VLS) Encoding
```
+----------+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
| SSSSSSSS |         STRING_DATA         |
+----------+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
* S: Size of string.
* STRING_DATA: String payload stored in S bytes.
```


## Header Format
```
+-----------------------------------------------------------------------------------+
|             "m3vcf" string in binary + 4 bytes for version (1.0.0.0)              |
+-----------------------------------------------------------------------------------+
| 01101101 00110011 01110110 01100011 01100110  XXXXXXXX XXXXXXXX XXXXXXXX XXXXXXXX |
+-----------------------------------------------------------------------------------+


+----------+-------------------------------------+VVVVVVVVVVVVVVVVVVVVVVVVVVVVV+
| PPPPPPPP | SSSSSSSS SSSSSSSS SSSSSSSS SSSSSSSS |     SAMPLE_ID_ARRAY ...     |
+----------+-------------------------------------+VVVVVVVVVVVVVVVVVVVVVVVVVVVVV+
* P: Ploidy level stored in 1 byte.
* S: Number of samples stored has VLQ.
* SAMPLE_ID_ARRAY: Array of length S that stores sample ID's in VLS encoding.
```


## Block Format
```
+----------+---------------------------------------+VVVVVVVVVVVVVVVVVVVV+
| MMMMMMMM | NNNNNNNN  NNNNNNNN  NNNNNNNN NNNNNNNN | SAMPLE_MAPPING ... |
+----------+---------------------------------------+VVVVVVVVVVVVVVVVVVVV+
+vvvvvvvvv+-------------------------------------------------------------------------+vvvvvvvvvv+vvvvvvvvv+vvvvvvvvv+VVVVVVVVVVVVVVVVVVVV+
|  CHROM  | PPPPPPPP PPPPPPPP PPPPPPPP PPPPPPPP PPPPPPPP PPPPPPPP PPPPPPPP PPPPPPPP |  MARKID  |   REF   |   ALT   | UNIQUE_HAP_ROW ... |
+vvvvvvvvv+-------------------------------------------------------------------------+vvvvvvvvvv+vvvvvvvvv+vvvvvvvvv+VVVVVVVVVVVVVVVVVVVV+
...


* M: Number of rows in unique haplotype matrix (number of markers in block).
* N: Number of columns in unique haplotype matrix.
* SAMPLE_MAPPING: Array of integers with a byte width of ceil(log2(N + 1) / 8). Array contains sample count times ploidy level elements.
* CHROM: Chromosome string stored has VLS.
* P: Chromosome pos stored in 8 bytes.
* MARKID: Marker ID string stored has VLS.
* REF: Reference haplotype stored has VLS.
* ALT: Alternate haplotype stored has VLS.
* UNIQUE_HAP_ROW: Array of '.', '0' or '1' bytes of size N.

```