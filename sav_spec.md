# Sparse Allele Vectors Specification

## Variable Length Integer (VLI) Encoding
All quantities are encoded in LEB128 format (https://en.wikipedia.org/wiki/LEB128). Encoded integers can start in the middle of a byte allowing 1 to 7 bits of data to prefix the integer.

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
+----------------------------------------------------------------+
|   "sav" string in binary + 4 bytes for version (Major.minor)   |
+----------------------------------------------------------------+
| 01110011 01100001 01110110 MMMMMMMM MMMMMMMM mmmmmmmm mmmmmmmm |
+----------------------------------------------------------------+
+vvvvvvvvv+~~~~~~~~~~~~~~+~~~~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVVVVVVVVVV+
|  CHROM  | PLOIDY_LEVEL | SAMPLE_SIZE |     SAMPLE_ID_ARRAY ...     |
+vvvvvvvvv+~~~~~~~~~~~~~~+~~~~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVVVVVVVVVV+
* CHROM: Chromosome string stored has VLS.
* PLOIDY_LEVEL: Ploidy level stored has VLI.
* SAMPLE_SIZE: Number of samples stored has VLI.
* SAMPLE_ID_ARRAY: Array of length SAMPLE_SIZE that stores sample ID's in VLS encoding.

+~~~~~~~~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVVVVVVVVVV+
| META_FIELDS_CNT |     META_FIELDS_ARRAY ...   |
+~~~~~~~~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVVVVVVVVVV+
* META_FIELDS_CNT: Number of metadata fields.
* META_FIELDS_ARRAY: Array of length META_FIELDS_CNT that stores metadata fiels in VLS encoding.

```

## Haplotype Pairs
Haplotype pairs are encoded in 1 or more bytes. The first bit of the first byte determines whether the pair represents a missing or alternate allele. The next 7 bits and any additional bytes in the pair make up a VLI that represents a zero-based offset from the previous non-zero haplotype.

### Example
```
An allele with an offset of 25.
+-+--------+
|0|001 1001|
+-+--------+

A missing haplotype with an offset of 8000.
+-+--------+ +---------+
|1|100 0000| |0111 1101|
+-+--------+ +---------+
```

## Record Format
```
+~~~~~~~~~+vvvvvvvvv+vvvvvvvvv+VVVVVVVVVVVVVVVVVVVVVV+~~~~~~~~~+VVVVVVVVVVVVVVVVVVVV+
|   POS   |   REF   |   ALT   | META_VALUE_ARRAY ... | HPA_SZ  | HAP_PAIR_ARRAY ... |
+~~~~~~~~~+vvvvvvvvv+vvvvvvvvv+VVVVVVVVVVVVVVVVVVVVVV+~~~~~~~~~+VVVVVVVVVVVVVVVVVVVV+

* POS: Chromosome pos stored has VLI.
* MARKID: Marker ID string stored has VLS.
* REF: Reference haplotype stored has VLS.
* ALT: Alternate haplotype stored has VLS.
* META_VALUE_ARRAY: Array of size META_FIELDS_CNT that stores metadata for each marker. Values correspond to META_FIELDS_ARRAY in header.
* HPA_SZ: Size of haplotype pair array stored as VLI with one bit prefix.
* HAP_PAIR_ARRAY: Array of size HPA_SZ that stores alleles or missing haplotypes in Haplotype Pair encoding.

```