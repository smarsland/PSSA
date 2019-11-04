# PSSA

Code to support the paper `Principal Symmetric Space Analysis' by Charles Curry, Stephen Marsland, and Robert McLachlan (Journal of Computational Dynamics 6(2), 2019).

riemanns2.m	 
  Computes the great circle that best approximates (in the Riemannian metric) data in S^2, together with the Karcher mean.
  Used to create Figure 2 of the paper.

spgo.m	
  Computes the best circle that best approximates (in the chordal metric) data in S^2.
  Used to create Figure 3 of the paper.
  
tori.m		
  Find the geodesic that best fits data in T^2.
  Used to create Figure 5 of the paper.
  
T13.m	
  Computes the best 1-torus and 2-torus for data in T^3 with fixed resonance relations.
  Used to create Figure 6 of the paper.

  
rotation_dataset.m	
  Computes the PSSA for data on S^2 x S^2 and plots the left branch
  Used to create Figure 7 of the paper.

torus_dataset.m
  Generates data on S^2 x S^2, finds the PSSA tree and plots the right branch.
  Used to create Figure 8 of the paper.

===

Helper Functions:
greatcircle.m		
  plots the great circle on S^2
frotgo.m and frotgo2.m		
  computes the PSSA tree for data (x,y) in S2 x S2. Draw the left and right branches respectively.
trueplot.m
  plots 3*N data matrix, where each column is a datapoint on the 2-sphere
drawsphere.m		
  draws the 2-sphere
drawtorus.m		
  draws the 2-torus
bestrot.m		
  finds the best rotation matrix for the left branch

