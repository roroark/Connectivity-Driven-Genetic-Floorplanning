## Connectivity Driven Genetic Floorplanning With Hard Macro Constraints

This repository contains code to automatically perform the [floorplanning](https://en.wikipedia.org/wiki/Floorplan_(microelectronics)) step of VLSI physical design, where we have to assign a shape and position to mutiple partitions (or sub-designs) having some soft-area which fit in any shape, and possibly  multiple hard macros (which fixed shape) constraining minimum dimensions. Partitions can interact with some paritions more than others (by have more pins/connections on an interface).

<img src="https://user-images.githubusercontent.com/18059416/91009876-1cc03080-e5ff-11ea-8fcd-497f843386a2.png" width="25%">

In general, floorplanning is an NP problem and requires exhaustive search to find the best solution. Genetic algorithms are one way to get a 'good enough' solution fast. Polish notation is used to encode floorplan information, however, which restricts us to only explore 'slicing floorplan' solutions. Floorplans where partitions that heavily interact with each other are placed beside each other are considered fitter and floorplans where partitions violating hard-macro requirements.

## Compile

Download all files into a folder, navigate to this folder and then compile "efp.cpp".
All dependencies are inluded.

## Usage

/path/to/efp.o <PARTITION_AREA_FILE> <PARTITION_CONNECTIVITY_FILE> <OUTPUT_FILE>
 
### Format of PARTITION_AREA_FILE
 
<NUMBER_OF_PARTITIONS>
for each of these partitions:
    <PARTITION_NAME> <SOFT_AREA> <HARD_MACRO_COUNT>
    for each of these hard macros:
        <HARDMACRO_WIDTH> <HARD_MACRO_HEIGHT>

For example:
3
BSS 122241 2
345 345
120 100
MSS 1234 0
TOP 421231 0

### Format of PARTITION_CONNECTIVITY_FILE:

NUMBER_OF_CONNECTIVITY_FIELDS
for each pair of partitions with non-zero interaction:
    <FROM_NAME> <TO_NAME> <CONNECTION_COUNT>

For example:
3
BSS TOP 200
MSS BSS 500
TOP MSS 420
