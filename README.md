## Theory :

Fourier-based field estimation is a Matlab toolbox to estimate magnetic field inhomogeneity due to spatial variation of magnetic susceptibility. 

In MRI it is assumed that the $B_0$ field is uniform, but due to differences in susceptibility in the body, the field will be non-uniform which creates artefacts. To compensate for the effect of the inhomogeneities in the field, shimming coils are installed in the scanner. To achieve the best possible correction with these coils, an accurate map of the magnetic field has to be obtained. That is exactly what this simulation does.

The derivation of the magnetic field equation is based on a Fourier transform of the susceptibility distribution and a Fourier transform of the point-dipole field. Only the z-component of the magnetic field $B_{dz}$ is important in the case of an external field $B_0$ in the z-direction. The derivation of the following equation for the FT of the magnetic field $\tilde B_{dz}$ is explained in ([Marques et al.](https://onlinelibrary.wiley.com/doi/10.1002/cmr.b.20034)).

$$ \tilde B_{dz} (\mathbf{k}) = \tilde M_{z} (\mathbf{k}) \cdot \mu_0 \bigg (\frac{1}{3} - \frac{k_z^2}{|\mathbf{k}|^2} \bigg) $$

With:

$$ |k|^2 = k_x^2 + k_y^2 + k_z^2 $$

If it is assumed that the susceptibility $\chi << 1 $, then an approximation for $M_{z} (\mathbf{r})$ can be made:

$$ M_{z} (\mathbf{r}) = \chi(\mathbf{r}) \frac{B_0}{\mu_0 (1 + \chi(\mathbf{r}))} \approx \chi(\mathbf{r}) \frac{B_0}{\mu_0} $$

The first equation can then be rewritten as:

$$ \tilde B_{dz} (\mathbf{k}) = \tilde \chi (\mathbf{k}) \cdot B_0 \bigg (\frac{1}{3} - \frac{k_z^2}{|\mathbf{k}|^2} \bigg) $$

For the zero frequency, $k=0$, this equation is undefined. The value for the zero frequency should be equal to the average field, for which some assumptions have to be made. In this toolbox the assumption is made that we are dealing with one of the following situations:
- a sphere with radius $a$ and susceptibility $\chi_i$ in an infinite medium of susceptibility $\chi_e$
- an infinitely long cylinder with the main axis parallel to $B_0$, radius $a$ and susceptibility $\chi_i$ in an infinite medium of susceptibility $\chi_e$

### Sphere in an infinite medium:

<p align="center">
<img src="https://user-images.githubusercontent.com/112189990/194596500-c4b6450d-8d6e-41f8-a768-fbed345f261e.png" width="200" height="230">
</p>
The derivations for the internal and external field of a sphere in an infinite medium are explained in (Brown et al.). This takes into account the Lorentz sphere correction, which accounts for the shift in the Larmor spin frequency inside the spherical body, when the local magnetic field is imaged. If the background material has a susceptibility of $\chi_e$, so not in vacuum, the BMS (background medium susceptibility) also has to be taken into account. This adds the term $\frac{1}{3} \chi_e B_0$, which also includes a Lorentz effect. 
<br/>
<br/>

- Internal field: $\frac{1}{3} \chi_e B_0$
- External field: $\frac{1}{3} (\chi_i - \chi_e) \cdot \frac{a^3}{r^3} (3 \cos^2(\theta) - 1) \cdot B_0 + \frac{1}{3} \chi_e B_0$

From this the average field value can be derived. For $r >> a$ , we can see that both the internal and external field will go to a value of $\frac{1}{3} \chi_e B_0$. Accordingly, the value for $k=0$ is set to this average.

$$\tilde B_{dz} (k=0) = \frac{1}{3} \chi_e B_0$$


### Infinitely long cylinder (|| $B_0$) in an infinite medium
<p align="center">
<img src="https://user-images.githubusercontent.com/112189990/194596320-76b668d3-5dbd-42f7-881e-e43b82f3653c.png" width="200" height="230">
</p>

The same Lorentz corrections and BMS shift apply in case of a cylinder, the derivations can again be found in (Brown et al.). 
<br/>
<br/>

- Internal field: $\frac{1}{6} (\chi_i - \chi_e) \cdot (3\cos^2(\theta) - 1) B_0 + \frac{1}{3} \chi_e B_0$
- External field: $\frac{1}{2} (\chi_i - \chi_e) \cdot \frac{a^3}{r^3} (3\cos^2(\theta) - 1) B_0 + \frac{1}{3} \chi_e B_0$

If $r>>a$ and $\theta = \frac{\pi}{2}$ (parallel to $B_0$), then the internal and external field again go to a value of $\frac{1}{3} \chi_e B_0$. 
## Overview :
The **ChiDist** class has 4 subclasses (**SheppLogan**, **Spherical**, **Cylindrical** and **Zubal**) and can be used to simulate the susceptibility distribution arising from a sphere, cylinder, Shepp-Logan or [Zubal phantom](http://noodle.med.yale.edu/zubal/data.htm). 

**FBFest.m** can be used to calculate the magnetic field offset produced by susceptibility distribution subject to a uniform external magnetic field B0 that is orientated along the z-axis.

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
