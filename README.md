# VMOMS
VMOMS Tokamak MHD Equilibrium Code

[VMOMS - A computer code for finding moment solutions to the Grad-Shafranov equation](https://doi.org/10.1016/0010-4655(82)90069-8)  
by L. Lao, R. M. Wieland, W. A. Houlberg and S. P. Hirshman.

The code was obtained from the [CPC archive](https://data.mendeley.com/journal/00104655)
and made to compile on a standard Linux system with `gfortran`.
VMOMS is archived there under the ID [`ABSH`](https://data.mendeley.com/datasets/g2z3562tv4/1).

## Building

```bash
> make
```

creates the `vmoms` executable.

## Running

The input file is `fort.10`.
Run the code as `./vmoms`.
The output data is in `fort.11`.
Some additional data (I think the flux surface geometry and the profiles) is in the plot output in `fort.12`.
