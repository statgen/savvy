# CVCF Specification

## Variable Length Quantity (VLQ) Encoding
Encoding is described at https://tools.ietf.org/html/rfc7541#section-5.1

## Variable Length String (VLS) Encoding
```
+~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
|   SIZE   |         STRING_DATA         |
+~~~~~~~~~~+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
* SIZE: Size of string encoded as VLQ.
* STRING_DATA: String payload stored in SIZE bytes.
```


## Header Format
```
+-------------------------------------------------------------------------+
|        "cvcf" string in binary + 4 bytes for version (1.0.0.0)          |
+-------------------------------------------------------------------------+
| 01100011 01110110 01100011 01100110 XXXXXXXX XXXXXXXX XXXXXXXX XXXXXXXX |
+-------------------------------------------------------------------------+
+~~~~~~~~~~~~~~+~~~~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVVVVVVVVVV+
| PLOIDY_LEVEL | SAMPLE_SIZE |     SAMPLE_ID_ARRAY ...     |
+~~~~~~~~~~~~~~+~~~~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVVVVVVVVVVV+
* PLOIDY_LEVEL: Ploidy level stored has VLQ.
* SAMPLE_SIZE: Number of samples stored has VLQ.
* SAMPLE_ID_ARRAY: Array of length SAMPLE_SIZE that stores sample ID's in VLS encoding.

```

## Haplotype Pairs
Haplotype pairs are encoded in 1 or more bytes. The first bit of the first byte is determines whether the pair represents a missing or alternate allele. The next 7 bits and any additional bytes in the pair make up a VLQ that represents an offset from the previous non-zero haplotype.

### Example
```
An allele with an offset of 100.
+-+-------+
|1|1100100|
+-+-------+

A missing haplotype with an offset of 2.
+-+-------+
|0|0000010|
+-+-------+
```

## Record Format
```
+vvvvvvvvv+~~~~~~~~~+vvvvvvvvvv+vvvvvvvvv+vvvvvvvvv+~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVV+
|  CHROM  |   POS   |  MARKID  |   REF   |   ALT   |  HPA_SZ  | HAP_PAIR_ARRAY ... |
+vvvvvvvvv+~~~~~~~~~+vvvvvvvvvv+vvvvvvvvv+vvvvvvvvv+~~~~~~~~~~+VVVVVVVVVVVVVVVVVVVV+


* CHROM: Chromosome string stored has VLS.
* POS: Chromosome pos stored has VLQ.
* MARKID: Marker ID string stored has VLS.
* REF: Reference haplotype stored has VLS.
* ALT: Alternate haplotype stored has VLS.
* HPA_SZ: Size of paplotype pair array stored as VLQ.
* HAP_PAIR_ARRAY: Array of size HPA_SZ that stores alleles or missing haplotypes in Haplotype Pair encoding.

```__