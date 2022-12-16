
## Theory 

The "Fourier-based field estimation" code allows one to estimate the magnetic field perturbation that arises when an object is placed within a magnetic field.

When an object is placed within an MRI scanner, it is assumed that the magnetic field experienced by the object is uniform and equal to the applied ( $B_0$ ) field. However, magnetic susceptibility ( $\chi$ ) differences between tissues in the body will lead to a non-uniform magnetic field within the body. This inhomogeneous magnetic field will cause image artefacts. These artefacts can be corrected for, if the magnetic field distribution is known. For this, an accurate map of the magnetic field must be acquired. This code allows the user to simulate magnetic fields, which can be useful for validating acquired field maps. 

In MRI, the $B_0$ field is aligned along the z-axis. When an object is placed within this field, it will become magnetized and only the z-component of the induced magnetization will be significant. The Fourier transform of z-component of the induced magnetic field can be expressed as follows (see [Marques et al.](https://onlinelibrary.wiley.com/doi/10.1002/cmr.b.20034) for a full derivation of this expression):

$$ \tilde B_{dz} (\mathbf{k}) = \tilde M_{z} (\mathbf{k}) \cdot \mu_0 \bigg (\frac{1}{3} - \frac{k_z^2}{|\mathbf{k}|^2} \bigg) $$ 

where the spatial frequency, $k$ is equal to $|k|^2=k_x^2+k_y^2+k_z^2$, $\mu_0$ is the permeability of free space, and $M_z$ is the induced magnetization along the z-axis and equal to:

$$ M_{z} (\mathbf{r}) = \chi(\mathbf{r}) \frac{B_0}{\mu_0 (1 + \chi(\mathbf{r}))} $$ 

If $\chi << 1 $, then we can approximate $M_{z} (\mathbf{r})$ as:

$$ M_{z} (\mathbf{r}) \approx \chi(\mathbf{r}) \frac{B_0}{\mu_0} $$

The first equation can then be rewritten as:

$$ \tilde B_{dz} (\mathbf{k}) = \tilde \chi (\mathbf{k}) \cdot B_0 \bigg (\frac{1}{3} - \frac{k_z^2}{|\mathbf{k}|^2} \bigg) $$

This equation allows us to simulate the magnetic field perturbation arising from a susceptibility distribution $\chi(r)$ when introduced within $B_0$. 

It should be noted that when $k=0$, the equation is undefined. $k=0$ is the spatial frequency with wavelength equal to zero, and $\tilde B_{dz} (\mathbf{k = 0})$ is otherwise interpreted as the average field. In order to avoid a singularity, one must assign a value to $\tilde B_{dz} (\mathbf{k} = 0)$, and for this, some assumptions must be made. 

### Setting the value of $\tilde B_{dz} (\mathbf{k} = 0)$ when the average magnetic field does not equal zero

In order to determine the appropriate value to assign to  $\tilde B_{dz} (\mathbf{k} = 0)$ we can consider two scenarios. 

#### Scenario 1: Sphere in an infinite medium

<p align="center">
<img src="https://user-images.githubusercontent.com/112189990/194596500-c4b6450d-8d6e-41f8-a768-fbed345f261e.png" width="200" height="230">
</p>

The derivation for the analytical solution of the magnetic field arising from a sphere placed within an infinite medium is given in Brown et al. This solution includes the Lorentz sphere correction. If the background material has a susceptibility of $\chi_e$ and sphere has a susceptibility of $\chi_i$, the magnetic field inside and outside of the sphere is expressed as:

- Internal field: $\frac{1}{3} \chi_e B_0$
- External field: $\frac{1}{3} (\chi_i - \chi_e) \cdot \frac{a^3}{r^3} (3 \cos^2(\theta) - 1) \cdot B_0 + \frac{1}{3} \chi_e B_0$

From this the average field value can be derived. For $r >> a$ , we can see that both the internal and external field will go to a value of $\frac{1}{3} \chi_e B_0$. $\tilde B_{dz} (\mathbf{k} = 0)$ can be set to $\frac{1}{3} \chi_e B_0$.

#### Scenario 2: Infinitely long cylinder in an infinite medium
<p align="center">
<img src="https://user-images.githubusercontent.com/112189990/194596320-76b668d3-5dbd-42f7-881e-e43b82f3653c.png" width="200" height="230">
</p>

The derivation for the analytical solution of the magnetic field arising from an infinite cylinder placed within an infinite medium is given in Brown et al. This solution includes the Lorentz sphere correction. If the background material has a susceptibility of $\chi_e$ and cylinder has a susceptibility of $\chi_i$, the magnetic field inside and outside of the cylinder is expressed as:

- Internal field: $\frac{1}{6} (\chi_i - \chi_e) \cdot (3\cos^2(\theta) - 1) B_0 + \frac{1}{3} \chi_e B_0$
- External field: $\frac{1}{2} (\chi_i - \chi_e) \cdot \frac{a^3}{r^3} \sin^2(\theta) \cos(2\phi) B_0 + \frac{1}{3} \chi_e B_0$

where $\theta$ is the angle between the direction of the main magnetic field and the central axis of the cylinder.

If $r>>a$, then the external field again goes to a value of $\frac{1}{3} \chi_e B_0$. Based on this, we can assume that  $\tilde B_{dz} (\mathbf{k} = 0) = \frac{1}{3} \chi_e B_0$.

### Setting the value of $\tilde B_{dz} (\mathbf{k} = 0)$ when the average magnetic field is equal to zero (i.e., a "demodulated" field)

Signals arising from an MRI scanner will be "demodulated". A consequence of this is that the average magnetic field within a measured field map is set to zero (here we call this a demodulated field) and any deviation from zero is due to susceptibility differences. 

<p align="center">
<img src="https://user-images.githubusercontent.com/112189990/206759060-6093c10d-b072-41ee-beb1-2eae9d184932.png" width="400" height="200">
</p>

In order to simulate this scenario, we can assume that $\tilde B_{dz} (\mathbf{k} = 0) = 0$. If the susceptibility differences between materials is known, then the demodulated field ( $\tilde B_{dz-demod} (\mathbf{k})$ ) can be computed as follows:

$$ \tilde B_{dz-demod} (\mathbf{k}) =  \Delta \tilde \chi (\mathbf{k}) \cdot B_0 \bigg (\frac{1}{3} - \frac{k_z^2}{|\mathbf{k}|^2} \bigg) $$


These final equations are the ones used in **FBFest**, which calculates the magnetic field offset produced by a susceptibility distribution subject to a uniform external magnetic field $B_0$ (oriented along the z-axis).

## Usage 

### Test script
Run the test script from the main folder (the folder containing FBFest), after adding it to the path.

A test script **test_calc_bdz** was developed for easy use of the FBFest function when testing with a spherical or cylindrical phantom. This test script allows a comparison to the analytical solutions for the sphere and cylinder, for which the equations are given in the theory. These equations are also adapted to give the solution for the demodulated field in ppm, so they only depend on the susceptibility difference and don't depend on the field strength of $B_0$. 

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
- **phi** (only for cylinder): angle between x-axis and the projection of the cylinder axis in the xy plane [rad], default values for phi are set and should not be changed:
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
Run the following commands from the main folder (the folder containing FBFest), after adding it to your path.

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
plot_along_axes(spherical_dBz_subsampled.volume, 'ppm')
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

