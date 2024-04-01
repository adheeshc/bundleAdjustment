# Bundle Adjustment

Using g2o on the BAL data set to demonstrate the BA experiments.
The BAL data set provides several scenes. The camera and landmark information in each scene are
given by a text file. This file stores the BA problem information in a line-by-line manner. For the detailed format,

see https://grail.cs.washington.edu/projects/bal.

We use the BALProblem class defined in common.h to read in the content of the file, and then use g2o to solve them.
