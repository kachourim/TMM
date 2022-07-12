# Transfer Matrix Method

Computes the power reflected and transmitted from a stack of dielectric or metallic layers using the transfer matrix method described in [1]. The implementation of this technique is based on the video series in [here](https://www.youtube.com/watch?v=l_CZFyJd4WI&list=PLLYQF5WvJdJVmCm4cDrKmek6cDJZWVomk).

The medium parameters may be dispersionless or dispersive. In that latter case, the parameters may be specified via a three-column txt file. The first column of the file is the wavelength (m), the second column is the real part of the relative permittivity/permeability (or the real part of the refractive index, *n*) and the third column is the imaginary part of the relative permittivity/permeability (or the imaginary part of the refractive index, *k*). See provided txt file (Au.txt) for illustration.

If no layer stack is specified, then the script simply computes the Fresnel coefficients due to scattering at the interface between the incidence and transmission media. In this case, the script does essentially the same thing as the one in [here](https://github.com/kachourim/FresnelCoefficients).

The different plotting options are described in the following sections, where *Rp* and *Tp* are the reflectance and transmittance for the parallel polarization (TM), and *Rs* and *Ts* are the reflectance and transmittance for the perpendicular polarization (TE). 

Note that this python script requires the `numpy`, `matplotlib` and `scipy` libraries. For the 2D computations below, the script is parallelized and uses all available cores (successfully tested on Debian 10 and Ubuntu 18.04). 

## 1D - Scattered power versus incidence angles at a single wavelength

To obtain the Fresnel coefficients for waves scattered at an interface between vacuum and a dielectric medium (assuming both media are dispersionless), you may use the following configuration:

```
# layer parameters
ER = []		# relative permittivity of each layer
MR = []		# relative permeability of each layer
L  = [] 	# thickness of each layer in m

# parameters of incidence medium
er_1 = 1
mr_1 = 1

# parameters of transmission medium
er_2 = 4
mr_2 = 1

# incidence angles in ° (theta can be a vector)
theta = linspace(0,90,500)
phi   = 0   

# free-space wavelength in m (can be a vector)
lam0 = 500e-9;
```

Note that the parameters `ER, MR and L` are empty (no stack) and that the exact value of the wavelength is not important in this case since the medium parameters are dispersionless. The result of this configuration is shown in the figure below. 

<img src="/images/1D_angle_noLayer.png" width="500">




## 1D - Scattered power versus wavelength at a single incidence angle

The script may also be used by specifying only one incidence angle and a set of wavelengths. Here is an example where the scattered power is computed at the interface between vacuum and gold for an incidence angle of 45°:

```
# layer parameters
ER = []	   	# relative permittivity of each layer
MR = []		# relative permeability of each layer
L  = [] 		# thickness of each layer in m

# parameters of incidence medium
er_1 = 1
mr_1 = 1

# parameters of transmission medium
er_2 = "Au.txt"
mr_2 = 1

# incidence angles in ° (theta can be a vector)
theta = 45
phi   = 0   

# free-space wavelength in m (can be a vector)
lam0  = linspace(200,900,500)*1e-9 

# Data file provides the permittivity/permeability (NK = 0) or the refractive index n/k (NK = 1)
NK = 1
```

Since gold is strongly dispersive within that wavelength range, the txt file containing the corresponding data is specified for `er_2`. Note that the script automatically interpolates the wavelength points not provided in the txt file. The resulting reflectance and transmittance spectra are plotted below.

<img src="/images/1D_wavelength.png" width="500">

## 1D - Anti-reflection coating

Here is an example of a single layer playing the role of an anti-reflection coating between vacuum and a dense dielectric medium. The refractive index of the layer is `n=sqrt(n1*n2)=2` and its thickness is `L = lambda/4`, for an operation wavelength of 500 nm. The corresponding configuration is

```
# layer parameters
ER = [2.0]	   			# relative permittivity of each layer
MR = [1.0]				# relative permeability of each layer
L  = [500e-9/(4*sqrt(2))] 		# thickness of each layer in m

# parameters of incidence medium
er_1 = 1
mr_1 = 1

# parameters of transmission medium
er_2 = 4
mr_2 = 1

# incidence angles in ° (theta can be a vector)
theta = 0
phi   = 0   

# free-space wavelength in m (can be a vector)
lam0  = linspace(50,2000,1000)*1e-9 
```

And the result is shown below.

<img src="/images/1D_coating.png" width="500">

## 2D - Surface plasmon

Example showing the dispersion of a surface plasmon for a  40 nm gold slab sandwiched between glass (incidence medium) and vacuum. 

```
# layer parameters
ER = ["Au.txt"]		# relative permittivity of each layer
MR = [1.0]			# relative permeability of each layer
L  = [40e-9] 			# thickness of each layer in m

# parameters of incidence medium
er_1 = 1.46**2
mr_1 = 1

# parameters of transmission medium
er_2 = 1
mr_2 = 1

# incidence angles in ° (theta can be a vector)
theta = linspace(40,60,200)
phi   = 0   

# free-space wavelength in m (can be a vector)
lam0  = linspace(400,800,200)*1e-9 

# Data file provides the permittivity/permeability (NK = 0) or the refractive index n/k (NK = 1)
NK = 1
```

The result is shown below, where we easily identify the surface plasmon dispersion in the reflection coefficient of the p-polarized (TM) wave.

<img src="/images/2D_SPP.png" width="750">




## 2D - Bloch surface wave

Finally, here is an example of the dispersion of a Bloch surface wave on a stack of dielectric layers, 5 * 2 layers (SiNx and SiO2),  sandwiched between glass and vacuum. The configuration file is

```
# layer parameters
ER = [(1.94-1j*1e-3)**2, (1.468-1j*1e-3)**2] * 5		# relative permittivity of each layer
MR = [1.0, 1.0] * 5						# relative permeability of each layer
L  = [250e-9, 460e-9] * 5					# thickness of each layer in m

# parameters of incidence medium
er_1 = 1.46**2
mr_1 = 1

# parameters of transmission medium
er_2 = 1
mr_2 = 1

# incidence angles in ° (theta can be a vector)
theta = linspace(30,60,200)
phi   = 0   

# free-space wavelength in m (can be a vector)
lam0  = linspace(1300,1800,200)*1e-9 
```
And the result is shown below, where we identify the dispersion of the Bloch surface wave in the reflection coefficient of the s-polarized (TE) wave.

<img src="/images/2D_BSW.png" width="750">


## Reference
[1] Rumpf, Raymond C. "Improved formulation of scattering matrices for semi-analytical methods that is consistent with convention." Progress In Electromagnetics Research 35 (2011): 241-261.
