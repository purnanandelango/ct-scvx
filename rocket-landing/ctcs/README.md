# Notes

## Convergence recipe

- Relax constraint bounds. Loosen scaling bounds of `bet`. Pick `trig_scl` close to 5.
- Successively warm-start and with harsher constraint parameters to obtain solution on coarser grid.
- Reduce `trig_scl` gradually to 0.5 and tighten scaling bounds of STC `bet` from [0,0.1] to [0,0.001]. 
- Tightening the bounds for `bet` scaling is especially important for good resolution of STC behavior when trigger function is close to active.