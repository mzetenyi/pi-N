# pi-N
> Codes for calculating dilepton production in pion-nucleon collisions.

This repository contains codes for the calculation of the angular distribution of di-electron production in pion-nucleon 
collisions using the model of [1]. In addition to calculating differential cross sections and anysotropy
parameters, the code can be use as an event generator, and can give the polarization density matrix of the virtual
photon decaying to the electron-positron pair. 
The model contains contributions of baryon resonances in the s- and u-channels.
Vector Meson Dominance is used: dileptons are created via an intermediate rho meson.

[1] E. Speranza, M. Zétényi, B. Friman, *"Polarization and dilepton anisotropy in pion-nucleon 
collisions"*, Phys. Lett. **B764** (2017) 282.

## Getting started
### Download

If you are not using `git` and `GitHub`, just click on the green "Clone or download" button in the top-right corner
of the [project page](https://github.com/mzetenyi/pi-N) and in the popup window
click "Download ZIP". Save the file and uncompress it.

### Compile

You need a c++11 compliant compiler and `make`. (I use g++ 4.8.5, and Gnu Make 4.1) If you don't have a suitable `g++`
compiler in your $PATH, you should modify line 20 of the [Makefile](https://github.com/mzetenyi/pi-N/blob/master/Makefile)
like this:
```
CPP = <ful path of compiler>
```
In order to compile, navigate to the 
main directory of the project (called `pi-N-master` if you downloaded the ZIP file and `pi-N` if you used `git` to 
obtain the files), and issue the command:
```
make piNdilep
```
You have to repeat this once more, because for some reason the linker does not find the freshly created object files.
(You will also get some warnings which you can ignore.) If the compilation was succesful you should have an `obj`
subdirectory containing object files, plus a `bin` directory containing the file `piNdilep`

### First run

Try running the code by entering
```
./jobs
```
This will take a few seconds, and then you will have a `results` directory with an `events` and a `density_matrix`
subdirectory, each containing an output file.

Take a look at the [`jobs`](https://github.com/mzetenyi/pi-N/blob/master/jobs) file! It's a shell script that first
creates the directories for storing the results and then runs the code twice. The first run generates 100 events
with a certain parameter set specified on the command line, the second run tabulates the density matrix for different
virtual photon masses and scattering angles, for the same set of parameters. The code writes to standard output, which is
redirected to files by the `jobs` script.

## Parameters
### Syntax

All parameters/options to the code can be specified in two ways:

1. as a command line option in the form `parameter=<value>` (no spaces are allowed),
2. in a file that is loaded by the command line option `load[<path_of_file>]`.

As an example for the latter, look at the file 
[dat/model_params](https://github.com/mzetenyi/pi-N/blob/master/dat/model_params).
As you can see, some of the parameters are arranged in groups via
```
group {
  param1 = <value1>   // comment
  param2 = <value2>
  ...
}
```
(Here spaces are allowed and any text after `//` is comment.) 

Parameters in a file can be overridden on the command line (the overriding options must be given after the load option!). 
In particular, grouped parameters can be overridden like this:
```
piNdilep group.param1=<new_value1> ...
```
E.g. `N1520.mass=1.49` will result in a calculation where the mass of the N(1520) resonance is artificially set
to 1.49 GeV.

Options without a value are also possible. These are typically switches that turn on and off various channels or decide
what quantities are to be calculated.

### Parameters of the code
#### Model parameters

The physical parameters of the model are conatined in the file 
[dat/model_params](https://github.com/mzetenyi/pi-N/blob/master/dat/model_params). Most importantly, this file
contains one group of parameters for each particle that is included in the model. Parameters include the
particle's spin, isospin, mass, width, branching ratios and various coupling constants. For baryon resonances, 
we give the orbital agular momentum of the corresponding partial wave instead of parity. All physical quantities
are given in powers of GeV. 

#### Swithes regulating various channels

Contribution of each baryon resonance can be included in the calculation by specifying its name on the command line.
The options `negN1440` or `negD1600` will multiply the corresponding amplitudes by -1. (This only effects the
interference terms with other resonance contributions.)
Including `sch` or `uch` in the command line results in a calculation of only s-channel (or only u-channel) contributions
of all resonances. (The default is calculate both s- and u-channels.) E.g.
```
piNdilep N1440 N1535 negN1535 ...
```
will include u- and s-channell contributions of N(1440) and N(1535), with an extra - sign between amplitudes of the
two resonances,
```
piNdilep N1675 uch ...
```
will calculate the u-channel contribution of N(1675).

#### Switches deciding what to calculate

The following switches specify the target of the calculation, exactly one of them has to be given:

* `event_gen`: generate events with a dilepton in the final state, list the electron and positron threemomenta.
* `density_matrix`: calculate the virtual photon polarizzation density matrix for specific photon mass and
                     scattering angle
* `tab_density_matrix`: tabulate the density matrix, varying values of the photon mass and the scattering angle
* etc.

##### `event_gen`

The parameter `Nevent=<val>` is used to specify the number of events.
