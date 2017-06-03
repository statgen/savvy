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
+~~~~~~~~~~~+vvvvvvvv+~~~~~~~~~~~~~~+~~~~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVV+
| BIT_WIDTH | CHROM  | PLOIDY_LEVEL | SAMPLE_SIZE | SAMPLE_ID_ARRAY ... |
+~~~~~~~~~~~+vvvvvvvv+~~~~~~~~~~~~~~+~~~~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVV+
* BIT_WIDTH: Bit width of genotype values (0-7). 
* CHROM: Chromosome string stored has VLS.
* PLOIDY_LEVEL: Ploidy level stored has VLI.
* SAMPLE_SIZE: Number of samples stored has VLI.
* SAMPLE_ID_ARRAY: Array of length SAMPLE_SIZE that stores sample ID's in VLS encoding.

+~~~~~~~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVV+
| FILE_INFO_SIZE | FILE_INFO_ARRAY ... |
+~~~~~~~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVV+
* FILE_INFO_SIZE: Number of key value pairs of file level info.
* FILE_INFO_ARRAY: Array of FILE_INFO_SIZE key value pairs of file level info.

+~~~~~~~~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVVVV+
| META_FIELDS_CNT | META_FIELDS_ARRAY ... |
+~~~~~~~~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVVVV+
* META_FIELDS_CNT: Number of metadata fields.
* META_FIELDS_ARRAY: Array of length META_FIELDS_CNT that stores variant level annotation fields in VLS encoding.

```

## Allele Pairs
Allele pairs are encoded in 1 or more bytes. The first BIT_WIDTH bits of the first byte encodes the non-zero value of the allele. The the remaining bits and any additional bytes in the pair make up a VLI that represents a zero-based offset from the previous non-zero allele.

Allele values can be decoded with the following formula (VALUE + 1) / 2 ^ BIT_WIDTH. With 1-bit values, 0.5 is considered a missing binary genotype. 2-bit values and greater represent binned posterior genotype probabilities.

### Example
```
A 1-bit missing allele with an offset of 25.
+-+--------+
|0|001 1001|
+-+--------+

A 1-bit alternate allele with an offset of 8000.
+-+--------+ +---------+
|1|100 0000| |0111 1101|
+-+--------+ +---------+
```

## Record Format
```
+~~~~~~~~~+vvvvvvvvv+vvvvvvvvv+VVVVVVVVVVVVVVVVVVVVVV+~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVVVV+
|   POS   |   REF   |   ALT   | META_VALUE_ARRAY ... | APA_SZ  | ALLELE_PAIR_ARRAY ... |
+~~~~~~~~~+vvvvvvvvv+vvvvvvvvv+VVVVVVVVVVVVVVVVVVVVVV+~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVVVV+

* POS: Chromosome pos stored has VLI.
* MARKID: Marker ID string stored has VLS.
* REF: Reference haplotype stored has VLS.
* ALT: Alternate haplotype stored has VLS.
* META_VALUE_ARRAY: Array of size META_FIELDS_CNT that stores metadata for each marker. Values correspond to META_FIELDS_ARRAY in header.
* APA_SZ: Size of allele pair array stored as VLI with one bit prefix.
* ALLELE_PAIR_ARRAY: Array of size APA_SZ that stores alternate alleles Allele Pair encoding.

```