## Overlaps for circular gaussians

The overlap integrals for two circular gaussian only depends on distance -- in fact is another circular gaussian.

For overlaps for two 40 km gaussians, the diameter is 56.5685424949238, the normalization is 0.00027583665414413394, and a and c are 0.001733132812499999

----------------------------------------
gaussian | diameter | norm | a and c 
---------|---------|----------|-------------
40 km base | 40   |0.00055167 | 0.003466
70 km base | 70   |0.00018013 | 0.001131
40 over 40 |56.56 |0.00027583 | 0.001733
40 over 70 |80.62 |0.00013580 | 0.000853
----------------------------------------------

This means that one can construct the matrices for the backus-gilbert method fairly quickly.

The SMAP L2C footprints are slightly elliptical (47 x 39km)

gaussian | major | minor | norm | a | c
-------|-------|-------|-------|-------|-------
SMAP to 70 | 84.31 | 80.13 | 0.000130646 | 0.0007814 | 0.000863732

So to calculate the overlap, we just need the relative NS and EW distance in km for the various source footprints relative to the target location.

