## Connectivity Driven Genetic Floorplanning With Hard Macro Constraints

This repository contains code to automatically perform the [floorplanning](https://en.wikipedia.org/wiki/Floorplan_(microelectronics)) step of VLSI physical design, where we have to assign a shape and position to mutiple partitions (or sub-designs) having some soft-area which fit in any shape, and possibly  multiple hard macros (which fixed shape) constraining minimum dimensions. Partitions can interact with some paritions more than others (by have more pins/connections on an interface).

<img src="https://user-images.githubusercontent.com/18059416/91009876-1cc03080-e5ff-11ea-8fcd-497f843386a2.png" width="25%">

In general, floorplanning is an NP problem and requires exhaustive search to find the best solution. Genetic algorithms are one way to get a 'good enough' solution fast. Polish notation is used to encode floorplan information, however, which restricts us to only explore 'slicing floorplan' solutions. Floorplans where partitions that heavily interact with each other are placed beside each other are considered fitter and floorplans where partitions violating hard-macro requirements.

## Compile

Download repository into a folder, navigate to this folder and then compile "efp.cpp".<br>
All dependencies are inluded.

## Usage

/path/to/efp.o <PARTITION_AREA_FILE> <PARTITION_CONNECTIVITY_FILE> <OUTPUT_FILE>
 
### Format of <PARTITION_AREA_FILE>
 
<NUMBER_OF_PARTITIONS><br>
for each of these partitions:<br>
&nbsp;&nbsp;&nbsp;&nbsp;<PARTITION_NAME> <SOFT_AREA> <HARD_MACRO_COUNT><br>
&nbsp;&nbsp;&nbsp;&nbsp;for each of these hard macros:<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<HARDMACRO_WIDTH> <HARD_MACRO_HEIGHT><br>

For example:<br>
3<br>
BSS 122241 2<br>
345 345<br>
120 100<br>
MSS 1234 0<br>
TOP 421231 0<br>

### Format of <PARTITION_CONNECTIVITY_FILE>

NUMBER_OF_CONNECTIVITY_FIELDS<br>
for each pair of partitions with non-zero interaction:<br>
&nbsp;&nbsp;&nbsp;&nbsp;<FROM_NAME> <TO_NAME> <CONNECTION_COUNT><br>

For example:<br>
3<br>
BSS TOP 200<br>
MSS BSS 500<br>
TOP MSS 420<br>

### Format of <OUTPUT_FILE>

for each partition specified in <PARTITION_AREA_FILE><br>
&nbsp;&nbsp;&nbsp;&nbsp;x1 y1 x2 y2


