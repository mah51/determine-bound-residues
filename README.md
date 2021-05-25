A simple python script used to determine residues within an antibody CDR that are bound to an antigen for a given antigen-antibody structural complex; presented as a PDB file.

CDR's are extracted using [cothia numbering scheme](http://www.chemogenomix.com/chothia-antibody-numbering), it is known that binding of 2 residues on either side of the CDR can occur, thus these are included (var EXTRA_RESIDUES controls this). CDR residues within 4Å of the antigen are classed as bound (var CONTACT_DISTANCE controls this). 
