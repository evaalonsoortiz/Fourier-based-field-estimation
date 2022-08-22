## Overview :
The **ChiDist** class has 4 subclasses (**SheppLogan**, **Spherical**, **Cylindrical** and **Zubal**) and can be used to simulate the susceptibility distribution arising from a sphere, cylinder, Shepp-Logan or [Zubal phantom](http://noodle.med.yale.edu/zubal/data.htm). 

**FBFest.m** can be used to calculate the magnetic field offset produced by susceptibility distribution subject to a uniform external magnetic field B0 that is orientated along the z-axis.

The folder **utils** includes a function to calculate the minimal distance between the region of interest and the edges of a matrix and a function that down-samples a 3D matrix by a chosen factor.

## References :

J.P. MARQUES, R. BOWTELL Concepts in Magnetic Resonance Part B (Magnetic Resonance Engineering), Vol. 25B(1) 65-78 (2005)

ZUBAL, I.G., HARRELL, C.R, SMITH, E.O, RATTNER, Z., GINDI, G. and HOFFER, P.B., Computerized three-dimensional segmented human anatomy, Medical Physics, 21(2):299-302 (1994)

## Dependencies :
[phantom3d.m](https://www.mathworks.com/matlabcentral/fileexchange/9416-3d-shepp-logan-phantom)

## Example usage :

To generate a susceptibility distribution with the shape of cylinder:

```
cylindrical_sus_dist = Cylindrical( matrix, image_res, radius, theta, susceptibility);
cylindrical_dBz = FBFest( cylindrical_sus_dist.volume, image_res, matrix );
