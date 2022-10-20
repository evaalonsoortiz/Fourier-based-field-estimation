## Theory :

Fourier-based field estimation is a Matlab toolbox to estimate magnetic field inhomogeneity due to spatial variation of magnetic susceptibility. 

In MRI it is assumed that the $B_0$ field is uniform, but due to differences in susceptibility in the body, the field will be non-uniform which creates artefacts. To compensate for the effect of the inhomogeneities in the field, shimming coils are installed in the scanner. To achieve the best possible correction with these coils, an accurate map of the magnetic field has to be obtained. That is exactly what this simulation does.

The derivation of the magnetic field equation is based on a Fourier transform of the susceptibility distribution and a Fourier transform of the point-dipole field. Only the z-component of the magnetic field $B_{dz}$ is important in the case of an external field $B_0$ in the z-direction. The following equation is used in this toolbox, the derivation is explained in ([Marques et al.](https://onlinelibrary.wiley.com/doi/10.1002/cmr.b.20034)).

$ \tilde B_{dz} (\mbold{k}) = \tilde M{z} (\mbold{k}) \cdot \mu_0 \bigg (\frac{1}{3} - \frac{k_z^2}{\mbold{k}} \bigg)


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
cylindrical_sus_dist = Cylindrical( matrix, image_res, radius, theta, [sus_in sus_out]);
cylindrical_dBz = FBFest( 'spherical', cylindrical_sus_dist.volume, image_res, matrix, sus_out);
