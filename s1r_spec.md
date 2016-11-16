# Sort-tile-recursive One-dimensional R-tree Specification

## Header Format
```
+----------------------------------------------------------------+
|   "s1r" string in binary + 4 bytes for version (Major.minor)   |
+----------------------------------------------------------------+
| 01110011 00110001 01110010 MMMMMMMM MMMMMMMM mmmmmmmm mmmmmmmm |
+----------------------------------------------------------------+

+------------------------------------------------------------------------------------+~~~~~~~~~+
| RRRRRRBB | NNNNNNNN NNNNNNNN NNNNNNNN NNNNNNNN NNNNNNNN NNNNNNNN NNNNNNNN NNNNNNNN | PADDING |
+------------------------------------------------------------------------------------+~~~~~~~~~+
R: Reserved bits.
B: The size of each node in bytes is specified as 8^(B+2). Possible bucket sizes are 64, 512, 4,096 and 32,768.
N: Record count.
P: Padding to end of current block.
```

## Internal Nodes
Internal nodes contain at most (block size / 16) entries. Entries are made up of two big-endian 64-bit integers that. The first of which represents the start value of the interval while the second represents the length of the interval.

## Leaf Nodes
Internal nodes contain at most (block size / 32) entries. Entries are made up of four big-endian 64-bit integers that. The first of which represents the start value of the interval. The second represents the length of the interval. The third represents the offset of the record in the target file. The fourth represents the length of the record in the target file. In most cases, the unit for the record offset and length is a byte.

## Tree Structure
The index file begins with the header block followed by the root block. Then each level of the tree going downward is stored from left to right. All nodes must be full except for the last node at each level.

The number of nodes at each level and the positions of child nodes are calculated using using the record count and block size.

Take for example a tree that contains 1,000,000 records with a node size of 4096 bytes. The leaf nodes fit a maximum of 128 entries, while the internal nodes fit 256. 1,000,000 is divided by 128 to get 7,813 leaf nodes. 7,813 is divided by 256 to get 31 internal nodes at the next level up. Since 31 is less than 256, the following level contains the root node.    

## Generating the Tree
Trees are created by first sorting all entries by the midpoint of their intervals. The tree is then built from the bottom up with the intervals of each entry in the internal nodes representing the smallest range that contains the child node's intervals. 

##

![title](s1r_diagram.png)