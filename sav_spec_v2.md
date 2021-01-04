
# SAV 2.0 Specification

Version 2 of the SAV file format extends the BCF file format at http://samtools.github.io/hts-specs/VCFv4.3.pdf (starting at page 27). The differences between BCF2 and SAV2 are as follows:
* The magic string starts with SAV instead of BCF.
* A type byte of zero indicates a sparse vector.
* A type byte with the fourth bit set indicates that PBWT is enabled for the vector.
* The sample size is not redundantly stored in each record. The most significant bit of the space used to store sample size in BCF records in used to indicate a PBWT reset. The rest of the bits are reserved.
* Files are compressed with blocked zstd instead of blocked gzip.
* Files use [S1R indices](./s1r_spec.md) instead of CSI, which are appended to the end of the SAV file instead of stored as a separate file.
 

The following reference is adapted from http://samtools.github.io/hts-specs/BCFv2_qref.pdf. The differences from the BCF quick reference are bolded.  

## **SAV2** Quick Reference 

In BCF2 **and SAV2**, each key in the FILTER, INFO and FORMAT fields is required to be defined in the VCF header.
For each record, a key is stored as an integer which is the index of its first appearance in the header. ‘PASS’
is always indexed at 0, which is special cased as VCF does not require the presence of this word.
In BCF2 **and SAV2**, a typed value consists of a typing byte and the actual value with type mandated by the typing
byte. In the typing byte, the lowest **three** bits give the atomic type. **The fourth bit specifies PBWT status.**  If the number represented by the higher
four bits is smaller than 15, it is the size of the following vector; if the number equals 15, the following typed
integer is the array size. The highest 4 bits of a Flag type equals 0 and in this case, no assumptions can be
made about the lower 4 bits. The table below gives the atomic types and their missing values:


| Bit 0–**2**| C type    | Missing value        | Description                          |
|-----------:|:----------|:---------------------|:-------------------------------------|
|         1  |  int8_t   |  0x80                |  signed 8-bit integer                |
|         2  |  int16_t  |  0x8000              |  signed 16-bit integer               |
|         3  |  int32_t  |  0x80000000          |  signed 32-bit integer               |
|       **4**|**int64_t**|**0x8000000000000000**|**signed 64-bit integer**             |
|         5  |  float    |  0x7F800001          |  IEEE 32-bit floating pointer number |
|         7  |  char     |  ‘\0’                |  character                           |

A **SAV2** file is **zstd** compressed and all multi-byte value are little endian.

### **Sparse Vector Type**
**A type set to zero indicates a sparse vector. Sparse vectors store only non-zero values and the offsets of those non-zero values. The offsets are relative to the preceding offset (0-based). For example, if the first 3 values of a vector are non-zero, then the corresponding offsets for those 3 values would be zero, zero, and zero. Similarly, if only the second and third values are non-zero, then the corresponding offsets for those 2 values would be one and zero.** 

**Sparse vectors encode the vector size in the same way as any other type, but it is followed by an additional byte OOOOVVVV where the four 'O' bits specify the offset type (possible type: 1-4) and the 'V' bits specify the value type (possible types: 1-5). It is then followed by a typed integer that specifies the sparse (non-zero) size of the vector. The data payload is then organized as an array of relative offsets followed by an array of non-zero values.**

**For example, A sparse float vector with a size of 10 that contains two non-zero values would be:**
**0xA0 0x00 0x01 0x05 0x01 0x02**
* **0xA0 is the size of the vector (10)**
* **0x00 is the sparse type code**
* **0x01 is the offset type (int8_t)**
* **0x05 is the float type (float)**
* **0x01 0x02 is a typed intger specifying the sparse size (2)**  

### **GT Encoding TODO!!!!**
A genotype (GT) is encoded as an integer vector with each integer describing an allele and its phase
w.r.t. the previous allele. The first allele does not carry the phase information. In the vector, each integer is
organized as ‘(allele+1)<<1|phased’ where allele is set to -1 if the allele in GT is a dot ‘.’ (thus the higher
bits are all 0). The vector is padded with missing values if the GT having fewer ploidy.


|Field | Description | Type | Value|
|------|-------------|------|------|
|magic | BCF2 magic string | char[5] | **SAV**\2\1 |
|l_text | Length of the header text, including any NULL padding | uint32 t | |
|text | NULL-terminated plain VCF header text | char[l text] | |
|*List of VCF records (until the end of the BGZF section)*|
|l_shared| Data length from CHROM to the end of INFO | uint32_t | |
|l_indiv | Data length of FORMAT and individual genotype fields | uint32_t | |
|CHROM | Reference sequence ID | int32_t | |
|POS | 0-based leftmost coordinate | int32_t| |
|rlen | Length of reference sequence | int32_t | |
|QUAL | Variant quality; 0x7F800001 for a missing value | float | |
|n_allele_info | n_allele<<16\|n_info | uint32_t | |
|n_fmt_sample | n_fmt<<24\|n_sample | uint32_t | |
|ID| Variant identifier |typed str| |
|*List of alleles in the REF and ALT fields (n=n allele)*|
|allele| A reference or alternate allele | typed str | |
|FILTER| List of filters; filters are defined in the dictionary |typed vec| |
|*List of key-value pairs in the INFO field (n=n info)*|
|info key | Info key, defined in the dictionary |typed int| |
|info value| Value |typed val| |
|*List of FORMATs and sample information (n=n fmt)*|
|fmt key| Format key, defined in the dictionary | typed int| |
|fmt type| Typing byte of each individual value, possibly followed by a typed int for the vector length | uint8_t+ | |
|fmt value | Array of values. The information of each individual is concatenated in the vector. Every value is of the same fmt type. Variable-length vectors are padded with missing values; a string is stored as a vector of char. | (by fmt_type) | |