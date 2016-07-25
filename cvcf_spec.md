# CVCF Specification

## Variable Length Integer (VLI) Encoding
All quantities are encoded in LEB128 format (https://en.wikipedia.org/wiki/LEB128). Encoded integers can start in the middle of a byte allowing 1 to 7 bits of data to prefix the integer. This prefixing is currently only being used with haplotype pairs.

## Variable Length String (VLS) Encoding
```
+~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
|   SIZE   |         STRING_DATA         |
+~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
* SIZE: Size of string encoded as VLI.
* STRING_DATA: String payload stored in SIZE bytes.
```


## Header Format
```
+-------------------------------------------------------------------------+
|        "cvcf" string in binary + 4 bytes for version (1.0.0.0)          |
+-------------------------------------------------------------------------+
| 01100011 01110110 01100011 01100110 XXXXXXXX XXXXXXXX XXXXXXXX XXXXXXXX |
+-------------------------------------------------------------------------+
+vvvvvvvvv+~~~~~~~~~~~~~~+~~~~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVVVVVVVVVV+
|  CHROM  | PLOIDY_LEVEL | SAMPLE_SIZE |     SAMPLE_ID_ARRAY ...     |
+vvvvvvvvv+~~~~~~~~~~~~~~+~~~~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVVVVVVVVVV+
* CHROM: Chromosome string stored has VLS.
* PLOIDY_LEVEL: Ploidy level stored has VLI.
* SAMPLE_SIZE: Number of samples stored has VLI.
* SAMPLE_ID_ARRAY: Array of length SAMPLE_SIZE that stores sample ID's in VLS encoding.

```

## Haplotype Pairs
Haplotype pairs are encoded in 1 or more bytes. The first bit of the first byte is determines whether the pair represents a missing or alternate allele. The next 7 bits and any additional bytes in the pair make up a VLI that represents an offset from the previous non-zero haplotype.

### Example
```
An allele with an offset of 25.
+-+--------+
|1|001 1001|
+-+--------+

A missing haplotype with an offset of 8000.
+-+--------+ +---------+
|0|100 0000| |0111 1101|
+-+--------+ +---------+
```

## Record Format
```
+~~~~~~~~~+vvvvvvvvvv+vvvvvvvvv+vvvvvvvvv+~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVV+
|   POS   |  MARKID  |   REF   |   ALT   |  HPA_SZ  | HAP_PAIR_ARRAY ... |
+~~~~~~~~~+vvvvvvvvvv+vvvvvvvvv+vvvvvvvvv+~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVV+

* POS: Chromosome pos stored has VLI.
* MARKID: Marker ID string stored has VLS.
* REF: Reference haplotype stored has VLS.
* ALT: Alternate haplotype stored has VLS.
* HPA_SZ: Size of paplotype pair array stored as VLI.
* HAP_PAIR_ARRAY: Array of size HPA_SZ that stores alleles or missing haplotypes in Haplotype Pair encoding.

```