
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


### Infinitely long cylinder in an infinite medium
<p align="center">
<img src="https://user-images.githubusercontent.com/112189990/194596320-76b668d3-5dbd-42f7-881e-e43b82f3653c.png" width="200" height="230">
</p>

The same Lorentz corrections and BMS shift apply in case of a cylinder, the derivations can again be found in (Brown et al.). 
<br/>
<br/>

- Internal field: $\frac{1}{6} (\chi_i - \chi_e) \cdot (3\cos^2(\theta) - 1) B_0 + \frac{1}{3} \chi_e B_0$
- External field: $\frac{1}{2} (\chi_i - \chi_e) \cdot \frac{a^2}{r^2} \sin^2(\theta) \cos(2\phi) B_0 + \frac{1}{3} \chi_e B_0$

If $r>>a$ and $\theta = 0$ (cylinder axis perpendicular to $B_0$), then the internal and external field again go to a value of $\frac{1}{3} \chi_e B_0$. 

### Assumptions :

$$ \tilde B_{dz} (\mathbf{k}) = \tilde \chi (\mathbf{k}) \cdot B_0 \bigg (\frac{1}{3} - \frac{k_z^2}{|\mathbf{k}|^2} \bigg) $$

with 

$$ \tilde B_{dz} (\mathbf{k = 0}) = \frac{1}{3} \chi_e B_0 $$

However, we often only know the susceptibility difference and not the individual $\chi_i$ and $\chi_e$ values. The calculations of the field can then still be done, assuming that we only want to know the frequency demodulated field in ppm. This is the same as what we measure in an MRI scan, and can be recognized by the field going to 0 ppm far away from the object. With this assumption, the resulting field only depends on the susceptibility difference $\chi_i - \chi_e$. We then use the susceptibility difference map (right) instead of the susceptibility distribution map (left).

<p align="center">
<img src="https://user-images.githubusercontent.com/112189990/206759060-6093c10d-b072-41ee-beb1-2eae9d184932.png" width="400" height="200">
</p>

Because the calculations are done in ppm (parts per million), the field also doesn't depend on the strength of the $B_0$ field. The equations for the frequency demodulated field ($\tilde B_{dz-demod} (\mathbf{k})) in ppm then reduce to the following: 

$$ \tilde B_{dz-demod}[ppm] (\mathbf{k}) = \tilde \Delta\chi (\mathbf{k}) \bigg (\frac{1}{3} - \frac{k_z^2}{|\mathbf{k}|^2} \bigg) \cdot 1e6 $$

with 

$$ \tilde B_{dz-demod}[ppm] (\mathbf{k = 0}) = \frac{1}{3} \cdot 0 \cdot B_0 \cdot 1e6 = 0 $$

These final equations are the ones used in **FBFest**, which calculates the magnetic field offset produced by a susceptibility distribution subject to a uniform external magnetic field $B_0$ (oriented along the z-axis).

## Usage :

### Test script
Run the test script from the main folder (the folder containing FBFest), after adding the subfolders (utils, test_scripts, external) to the path.

A test script **test_calc_bdz** was developed for easy use of the FBFest function when testing with a spherical or cylindrical phantom. This test script allows a comparison to the analytical solutions for the sphere and cylinder, for which the equations are given in the theory. These equations are also adapted to give the solution for the frequency demodulated field in ppm, so they only depend on the susceptibility difference and don't depend on the field strength of $B_0$. 

Three flags in the beginning of the test script give the user some choices for the simulation. 
- **phantom**: the choice between "sphere" or "cylinder"
- **field**: the choice between "demodulated" or "offset". The default is demodulated. When only the susceptibility difference is known we are restricted to the default, however when $\chi_i$ and $\chi_e$ are known separately, the simulation of the field offset (not frequency demodulated) can be done. 
- **unit**: the choice between "ppm" or "Hz". Default is ppm. If ppm is chosen, then the simulation does not depend on the strength of $B_0$, if Hz is chosen then it does. 

The susceptibility distribution is created using the **ChiDist** class, which has 4 subclasses (**SheppLogan**, **Spherical**, **Cylindrical** and **Zubal**). In the test script, the Spherical and Cylindrical subclasses are used to simulate the field offset from a spherical or cylindrical phantom. The subclasses SheppLogan and Zubal can be used to simulate the field offset from the Shepp-Logan phantom or the [Zubal phantom](http://noodle.med.yale.edu/zubal/data.htm). For these phantoms there is no analytical solution. 

To create the susceptibility distributions, some parameters have to be set:
- **radius**: radius of the sphere or cylinder [mm]
- **theta** (only for cylinder): angle of rotation of the cylinder axis around y-axis [rad]
-   - theta = 0: cylinder axis parallel to z-axis and $B_0$
    - theta = $\pi/2$: cylinder axis perpendicular to z-axis and $B_0$
- **phi** (only for cylinder): angle between x-axis and measurement axis in the xy plane [rad], default values for phi are set and should not be changed:
    - phi_x = 0 (measurement along x-axis)
    - phi_y = $\pi/2$ (measurement along y-axis)
- **susin**: value of $\chi_i$ [ppm], only use this when you want to calculate the field offset
- **susout**: value of $\chi_e$ [ppm], only use this when you want to calculate the field offset
- **sus_diff**: value of $\chi_i - \chi_e$ [ppm], this value is used by default for the simulation of the demodulated and field offset. NB: if you have the specific values of $\chi_i$ and $\chi_e$ then you can change this parameter to susin - susout
- **dim_without_buffer**: field of view (FOV) in 3 dimensions (y, x, z)
- **dim**: FOV with buffer in 3 dimensions (y, x, z). The phantom is padded with zeros to the dimensions of dim, this improves the result from the Fourier Transform in FBFest, and will give a more accurate result of the simulated field. 
- **res**: resolution of the phantom [mm] in 3 dimensions, this refers to the physical size of a voxel in the FOV. Default is 1mm.
- **subsample_factor**: factor with which the simulated field is subsampled to match the resolution of the scanner. Default is 1 (no subsampling), a factor of 2 will give a resolution of 2mm if res is set to 1mm. The subsampling is done using the function **subsample_3D** in the utils folder.


The susceptibility distribution is then made using the ChiDist subclasses:
- **sphere**: sus_dist_volume = Spherical(dim_without_buffer, res, radius, [sus_diff 0]).volume;
- **cylinder**: sus_dist_volume = Cylindrical(dim_without_buffer, res, radius, theta, [sus_diff 0]).volume;
- **Zubal**: sus_dist_volume = Zubal('modified_zubal.nii').volume;
<br/>

### From command line
Run the following commands from the main folder (the folder containing FBFest), after adding the subfolders (utils, test_scripts, external) to the path.

#### 1. Spherical phantom
```
spherical_sus_dist = Spherical(matrix, image_res, radius, [sus_diff 0]);
spherical_dBz = FBFest('spherical', spherical_sus_dist.volume, image_res, matrix);
```

Example:

```
spherical_sus_dist = Spherical([128 128 128], [1 1 1], 10, [9 0]);
spherical_dBz = FBFest('spherical', spherical_sus_dist.volume, spherical_sus_dist.image_res, spherical_sus_dist.matrix);
```

Visualization:

```
plot_along_axes(spherical_dBz.volume, 'ppm')
sliceViewer(spherical_dBz.volume,'Colormap',colormap(parula),'SliceDirection',[1 0 0]); % or direction [0 1 0] or [0 0 1]
slice(spherical_dBz.volume, 64, 64, 64)
```

#### 2. Cylindrical phantom

```
cylindrical_sus_dist = Cylindrical(matrix, image_res, radius, theta, [sus_diff 0]);
cylindrical_dBz = FBFest('cylindrical', cylindrical_sus_dist.volume, image_res, matrix);
```

Example:

```
cylindrical_sus_dist = Cylindrical([128 128 128], [1 1 1], 10, pi/2, [9 0]);
cylindrical_dBz = FBFest('cylindrical', cylindrical_sus_dist.volume, cylindrical_sus_dist.image_res, cylindrical_sus_dist.matrix);
```

Visualization:

```
plot_along_axes(cylindrical_dBz.volume, 'ppm')
sliceViewer(cylindrical_dBz.volume,'Colormap',colormap(parula),'SliceDirection',[1 0 0]); % or direction [0 1 0] or [0 0 1]
slice(cylindrical_dBz.volume, 64, 64, 64)
```

#### 3. Zubal phantom
Important: first add the nifti file containing the modified Zubal phantom to the main folder
```
zubal_sus_dist = Zubal('zubal_modified.nii');
zubal_dBz = FBFest('Zubal', zubal_sus_dist.volume, zubal_sus_dist.image_res, zubal_sus_dist.matrix);
```

Visualization:

```
plot_along_axes(zubal_dBz.volume, 'ppm')
sliceViewer(zubal_dBz.volume,'Colormap',colormap(parula),'SliceDirection',[1 0 0]); % or direction [0 1 0] or [0 0 1]
slice(zubal_dBz.volume, 64, 64, 64)
```

#### 4. Buffer

```
spherical_sus_dist = Spherical(matrix, image_res, radius, [sus_diff 0]);
spherical_buffer_dBz = FBFest('spherical', spherical_sus_dist.volume, image_res, matrix, **buffer_dim**);
```

Example:

```
spherical_sus_dist = Spherical([128 128 128], [1 1 1], 10, [9 0]);
spherical_buffer_dBz = FBFest('spherical', spherical_sus_dist.volume, spherical_sus_dist.image_res, spherical_sus_dist.matrix, [256 256 256]);
```

Visualization:

```
plot_along_axes(spherical_buffer_dBz.volume, 'ppm')
```

#### 5. Subsample

```
spherical_sus_dist = Spherical(matrix, image_res, radius, [sus_diff 0]);
spherical_dBz = FBFest('spherical', spherical_sus_dist.volume, image_res, matrix);
spherical_dBz_subsampled = sub_sample_3D(spherical_dBz.volume, factor);
```

Example:

```
spherical_sus_dist = Spherical([128 128 128], [1 1 1], 10, [9 0]);
spherical_dBz = FBFest('spherical', spherical_sus_dist.volume, spherical_sus_dist.image_res, spherical_sus_dist.matrix);
spherical_dBz_subsampled = sub_sample_3D(spherical_dBz.volume, [2 2 2]);
```

Visualization:

```
plot_along_axes(spherical_dBz_subsampled, 'ppm')
```

#### 6. Get results in Hz instead of ppm

```
spherical_sus_dist = Spherical(matrix, image_res, radius, [sus_diff 0]);
spherical_dBz_ppm = FBFest('spherical', spherical_sus_dist.volume, image_res, matrix);
spherical_dBz_Hz = spherical_dBz_ppm.volume * B0[T] * Larmor_frequency[MHz/T];
```

Example:

```
spherical_sus_dist = Spherical([128 128 128], [1 1 1], 10, [9 0]);
spherical_dBz_ppm = FBFest('spherical', spherical_sus_dist.volume, spherical_sus_dist.image_res, spherical_sus_dist.matrix);
spherical_dBz_ppm_vol = spherical_dBz_ppm.volume;
spherical_dBz_Hz_vol = spherical_dBz_ppm_vol* 3 * 42.5775;
```


Visualization:

```
plot_along_axes(spherical_dBz_ppm_vol, 'ppm')
plot_along_axes(spherical_dBz_Hz_vol, 'Hz')
```


#### 7. Get field offset instead of frequency demodulated field

```
spherical_sus_dist = Spherical(matrix, image_res, radius, [sus_diff 0]);
spherical_dBz_demod = FBFest('spherical', spherical_sus_dist.volume, image_res, matrix);
spherical_dBz_offset = spherical_dBz_demod.volume + susout / 3;
```

Example:

```
spherical_sus_dist = Spherical([128 128 128], [1 1 1], 10, [9 0]);
spherical_dBz_demod = FBFest('spherical', spherical_sus_dist.volume, spherical_sus_dist.image_res, spherical_sus_dist.matrix);
spherical_dBz_demod_vol = spherical_dBz_demod.volume;
spherical_dBz_offset_vol = spherical_dBz_demod_vol + 1/3;
```


Visualization:

```
plot_along_axes(spherical_dBz_demod_vol, 'ppm')
plot_along_axes(spherical_dBz_offset_vol, 'ppm')
```


## References :

J.P. MARQUES, R. BOWTELL Concepts in Magnetic Resonance Part B (Magnetic Resonance Engineering), Vol. 25B(1) 65-78 (2005)

BROWN, W.B., CHENG, Y-C.N., HAACKE, E.M., THOMPSON, M.R. and VENKATESAN, R., Magnetic resonance imaging : physical principles and sequence design, chapter 25 Magnetic Properties of Tissues : Theory and Measurement. John Wiley & Sons, 2014.

ZUBAL, I.G., HARRELL, C.R, SMITH, E.O, RATTNER, Z., GINDI, G. and HOFFER, P.B., Computerized three-dimensional segmented human anatomy, Medical Physics, 21(2):299-302 (1994)

## Dependencies :
[phantom3d.m](https://www.mathworks.com/matlabcentral/fileexchange/9416-3d-shepp-logan-phantom)

