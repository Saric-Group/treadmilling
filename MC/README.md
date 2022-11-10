Description of the code in this folder:

Code to generate MC-MD simulations of treadmilling filaments. 
Here polymerisation, depolymerisation and nucleation are implemented as separate and independent processes.
Each is governed by a rate of occurrence K_{on/off/nuc} [steps^{-1}]
Initial conditions comprise a group of FtsZ filaments and a group of free individual monomers distributed on the membrane.
This is characterised by a list of heads, tails and free monomers (each with their growing direction).
For each "active" monomer (head, tail or free) we compute (probabilistically) the next time it will be involved in a treadmilling event:
	- Heads: from K_{on} we find t_{on}, the time it will polymerise a new monomer onto itself
	- Tails: from K_{off} we find t_{off}, the time it will depolymerise from its current filament
	- Free: from K_{nuc} we find t_{nuc}, the time it will nucleate a new filament by polymerising a new monmer onto itself
All this information (monomers ids, growth directions and event times) is stored in dynamic arrays for each type

NUCLEATION EVENT:
- From event time find id and dir of involved monomer (free)
- Choose polymerising monomer (nid) from available pool of free monomers
- Find the new position and check availability of the nucleation move
- if ACCEPTED:
	+ Update HEADS, TAILS & FREE arrays: 
		> add nid to HEADS with direction dir (from id) and a newly computed t_{on} (from K_{on})
		> add id to TAILS with its direction (dir) and a newly computed t_{off} (from K_{off})
		> remove id & nid from FREE
	+ Change position of nid
	+ Add a bond between id & nid in BONDS MATRIX
- Find next event: time & type
- if time for next event > current time:
	+ Write new input for LAMMPS
	+ Run until next event

POLYMERISATION EVENT:
- From event time find id and dir of involved monomer (heads)
- Choose polymerising monomer (nid) from available pool of free monomers
- Find the new position and check availability of the polymerisation move
- if ACCEPTED:
	+ Update HEADS & FREE arrays: 
		> replace id for nid in HEADS with direction dir (from id) and a newly computed t_{on} (from K_{on})
		> remove nid from FREE
	+ Change position of nid
	+ Add a bond between id & nid in BONDS MATRIX
- Find next event: time & type
- if time for next event > current time:
	+ Write new input for LAMMPS
	+ Run until next event

DEPOLYMERISATION EVENT:
- From event time find id and dir of involved monomer (tails)
- Find neighbour (from BONDS MATRIX) which will become the new tail: nid
- ALWAYS:
	+ Update HEADS & FREE arrays: 
		> replace id for nid in TAILS with direction dir (from id) and a newly computed t_{off} (from K_{off})
		> add id to FREE
	+ Remove bond between id & nid in BONDS MATRIX
- Find next event: time & type
- if time for next event > current time:
	+ Write new input for LAMMPS
	+ Run until next event

Flowchart:

 ---------------
|				|
|				|
|	   I.C.		|
|				|
|				|
 ---------------

		||
		\/

 ---------------
|				|
|				|
|	   	MD		|
|				|
|				|
 ---------------

		||
		\/

 ---------------			 -------------------			 ---------------
|				|			|					|			|				|
|				|			|	  MC event: 	|			|				|
|	STATE t0	|	--->	|  nucl / on / off  |	--->	|	STATE t1	|
|				|			|					|			|				|
|				|			|					|			|				|
 ---------------			 -------------------			 ---------------

		/\															||
		||															\/
		||
		||													 ---------------
		||													|				|
		||													|				|
		|| <----<-------<-------<-------<-------<-------<-- |	   	MD		|
															|				|
															|				|
															 ---------------
