The ChiDist class has 3 subclasses (SheppLogan, Spherical and Cylindrical) and can be used to simulate the susceptibility distribution arising from a sphere, cylinder or Shepp-Logan phantom. 

fourier_based_field_est.m can be used to calculate the magnetic field offset produced by susceptibility distribution subject to a uniform external magnetic field B0 that is orientated along the z-axis.

REFERENCES: J.P. MARQUES, R. BOWTELL Concepts in Magnetic Resonance Part B (Magnetic Resonance Engineering), Vol. 25B(1) 65?78 (2005)

Requirements: 
[phantom3d.m](https://www.mathworks.com/matlabcentral/fileexchange/9416-3d-shepp-logan-phantom)

Example usage:

```To generate a susceptibility distribution with the shape of cylinder:

my_sus_dist = Cylindrical( matrix, image_resolution, radius, theta, susceptibility);

my_sus_dist.save('cylindrical_sus_dist.nii');

Bdz = fourier_based_field_est('cylindrical90_sus_dist.nii');```
