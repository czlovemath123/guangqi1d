The monochromatic dust opacity \kappa_lambda [cm^2/g] computed by 
Semenov et al. 2003 for the case of the "iron-poor" silicate mineralogy. 
The folder names can be translated as follows:

Comp_aggregates - composite aggregate model of dust grains,
Hom_aggregates - homogeneous aggregate model of dust grains,
Hom_spheres - homogeneous spherical grains,
Comp_spheres - composite spherical grains,
Multishell_spheres - multishell spherical grains,
Porous_comp_spheres - porous composite spherical grains,
Porous_multishell_spheres - porous multishell spherical grains.

Each folder contains opacity data files for erey temperature region.
That is, "1.dat" corresponds to the first temperature region and so on.

Dear Zhuo,

Thank you for considering my dust opacities! Indeed, the files could have been more
clearly described. The 5 columns in e.g. the "1.dat" file for comp. aggregates are as follows:
Wavelength [microns], Kextm, Kext, Kabs, Ksca (all 3 in cm^2/g).
Here, Kextm=<Qabs+Qsca> * (1-<Qsca*<cos>>/<(Qabs+Qsca)>),
Kext =<Qext> = <Qabs> + <Qsca>,
Kabs =<Qabs>,
and Ksca =<Qsca>.

If you are mainly interested in (sub)-mm opacities, then you can use
Kabs or Kext (Kext almost equal to Kabs at long wavelengths, as scattering
and thus Ksca is negligible). If you need opacities at shorter wavelengths, you 
have to choose between total opacities (Kext), or absorption only
(Kabs) or scattering only (Ksca).

Regarding T-ranges, these boundaries are not that precisely defined,
and one has to take into account that the transitions are likely not
sharp. Thus, the difference between each temperature range is a transition
region where opacities drop because one of the key solid material begins
to evaporate (e.g., water ice at ~150 K). 

Hope that helps,

Dima
