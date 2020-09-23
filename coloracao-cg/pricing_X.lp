\ SCIP STATISTICS
\   Problem name     : pricing
\   Variables        : 3 (3 binary, 0 integer, 0 implicit integer, 0 continuous)
\   Constraints      : 2
\   Obj. scale       : 1
\   Obj. offset      : 0
Maximize
 Obj: +1 x_0 +1 x_1
Subject to
 edge0: +1 x_0 +1 x_1 <= +1
 edge1: +1 x_0 +1 x_2 <= +1
Bounds
 0 <= x_0 <= 1
 0 <= x_1 <= 1
 0 <= x_2 <= 1
Binaries
 x_0 x_1 x_2
End
