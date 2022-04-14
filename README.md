# bsavoie-MD_tool
in-group MD tool for you to do various MD analysis on lammps trajectory

This repo is originally written by: Aditi, Bumjoon, Dylan G, Jack, Lin (in alphabetical order)

This document is to serve as a context and guide for continuing the script development. The goal is/was to incorporate object oriented programming and have common scripts in one place which perform MD simulation analysis.                                                                                                                                                                               

**********************************************************
# CONTEXT
**********************************************************

As a starting point we introduced classes and used them in evaluating RDFs for a system. In that we wrote 5 scripts:
- mol_classes.py: Defines the classes.
- data_parser.py: Reads data file.
- frame_generator.py: Reads trajectory file.
- rdf.py: Parses the RDF
- rdf_parallel.py: Runs rdf.py parallely for efficiency.
You can refer to individual scripts for documentation on how to use the scripts. Additionally, to maintain consistent documentation, the documentation guidelines are written in FORMAT_GUIDELINES


**********************************************************
# FUTURE DEVELOPMENT
**********************************************************

In future, to continue, three things could be beneficial:
1. Integrating this with some of the functionality in TAFFI
2. Writing other analysis scripts using this
3. Developing the 5 scripts further (It can be found in the scripts itself)

**********************************************************
### INTEGRATING WITH TAFFI
**********************************************************
Example:
TAFFI can determine distinct molecule types for a given system. This is very useful for a mixture system, when you want to analyze properties of a specific molecule.
This functionality can be found in function parse_molecules in ~/bin/taffi/Parsers/msd_parse_multi.py. Given data file(s), it finds distinct molecules by identifying disconnected subgraphs using adjacency matrix.
Then assigns each distinct molecule a type based on the neighbors and their charges.
This list of molecule types can be added as an attribute to the AtomList class. This can be either done in data_file_parser.py or as a method within the AtomList class in mol_classes.py.

Useful TAFFI functionality which can be integrated:
1. Molecule type in function parse_molecules in *~/bin/taffi/Parsers/msd_parse_multi.py*
2. Adjacency matrix generator given coordinates— Table_generator function in *~/bin/taffi/Lib/adjacency.py*

**********************************************************
### Rewriting other analysis scripts using this   
**********************************************************
Example:
Parsing MSD for atoms of a given atom type "t". There is an existing script *~/bin/taffi/Parsers/msd_parse_multi.py.
The MSD parser follows the basic skeleton:
- Parsing data file-- use data_parser.py

- Identify the IDs of atoms which have the atom type "t"-- get_idx method in AtomList class

- Read trajectory using-- use frame_generator

- Calculate the MSD--  the MSD calculation part from ~/bin/taffi/Parsers/msd_parse_multi.py

Analysis script wishlist:
1. MSD parser
2. Autocorrelation parser
