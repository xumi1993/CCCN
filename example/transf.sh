#!/bin/bash
sac<<eof
r IC.SSE.00.BHZ
rmean;rtr
lp c 0.4
interp delta 1
transfer FROM POLEZERO SUBTYPE SAC_PZs_IC_SSE_BHZ_00 TO VEL freq 0.005 0.01 0.3 0.35
rmean;rtr
bp c 0.02 0.067
w over
q
eof
