simulate - Boilerplate Molecular Dynamics Code
==============================================

This is an experimental molecular dynamics simulation code developed to 
resemble the MOLDYN code developed at the Interfacial Molecular Science 
Laboratory at Rutgers University (IMSL; see http://imslab.rutgers.edu/).
Rewritten in C99, it was written with two goals in mind:

1. Function similarly to the original FORTRAN77 code from the user perspective 
   yet being an unencumbered, clean re-implementation that allows for the 
   integration of newer, high-performance features and APIs such as CUDA, 
   OpenACC, and automatic vectorization
2. Provide a freely available implementation of the Dissociative Water Potential
   (DWP) developed at the Interfacial Molecular Science Laboratory that 
   reproduces the behavior and properties reported by the original publications
   describing the model.

License
-------
This software is made available under the Creative Commons Attribution-
NonCommercial-ShareAlike license.  See the LICENSE file for more details.
