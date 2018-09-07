# EoL-Cloth
Eulerian-on-Lagrangian Cloth Simulation
## Info
EoL-Cloth is a simulation tool for simulating the phenomenon of cloth sliding over sharp geometric features such as boxes and points. It was developed at Texas A&M University and the coinciding publication was presented at SIGGRAPH 2018.  

Check out the project page [here](http://faculty.cs.tamu.edu/sueda/projects/eol-cloth/ "EoL-Cloth").

## Installing
#### Dependencies
* [CMake](https://cmake.org/ "CMake")
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page "Eigen")

Technically the code will run with just these two libraries, but simulations must be collision free and nothing will render.

* [GLM](https://glm.g-truc.net/0.9.9/index.html "GLM")
* [GLFW](http://www.glfw.org/ "GLFW")
* [GLEW](http://glew.sourceforge.net/ "GLEW")

To run the simulation "online" with rendered output these three libraries must be installed

* [MOSEK](https://www.mosek.com/ "Mosek")
OR
* [GUROBI](http://www.gurobi.com/ "Gurobi")

To handle collisions and in turn the EoL dynamics of the simulation, a quadratic program (QP) must be used. The Mosek and Gurobi libraries are currently supported using corresponding QP wrappers from [here](http://www.google.com/ "Mosek QP Wrapper") and [here](https://github.com/jrl-umi3218/eigen-gurobi "Gurobi QP Wrapper") respectively. These wrappers do not need to be installed seperately, they are already a part of the EoL-Cloth src.

While it does not need to be seperetely installed, we use a modified version of [ARCSim](http://graphics.berkeley.edu/resources/ARCSim/ "ARCSim") to handle our remeshing. Feel free to look into some of the great work that has been done with it. 

#### Building
```sh
git clone https://github.com/sueda/eol-cloth
cd eol-cloth
mkdir build
```
Before building with CMake, open the CMakeLists.txt file with a text editor (I know this is bad practice). While the build options can be set when we run CMake, for each option turned ON the corresponding dependency file path must be set. I do this for ease of use with multiple external library versions installed. Set the paths as shown to match the machine. Note some of the binaries are defined for the solvers which are platform specific. 
```sh
cd build
cmake [options] ..
make
```

Where the options, if not already changed in the CMakeLists.txt, are:
 * `-DCMAKE_BUIlD_TYPE=Release` Build in Release mode
 * `-DONLINE=ON` default is `OFF`
 * `-DMOSEK=ON` default is `OFF`
 * `-DGUROBI=ON` default is `OFF`

## Running
```sh
./eol-cloth <generalSettings.json> <simulationSettings.json>
```
Settings are defined using json files and their file paths are passed as the two command line arguments.
#### General Settings
* `online` : `true/false`
* `exportObjs` : `true/false`
* `RESOURCE_DIR` : path the the directories with rendering shaders
* `OUTPUT_DIR` : path to export location
 
Fairly self explanatory.
#### Simultion Settings
Most of these settings are described better in the actuall example simulationSettings json files, but as a quick overview:
* `solver` : `none/mosek/gurobi`
* `remeshing` : `true/false` which equates to on/off
* `EOL` : `true/false` which equates to on/off
* `Cloth` : These settings cover initial cloth `shape, resolution, and position`, `materials`, `remeshing parameters`, and `fixed points`
* `Obstacles` : These settings cover `collision threshold` and a basic definition structure for building `points` and `boxes`
 
## Exporting 
To export our objects we use an in lab developed tool we call Brender. After defining an export directory and turning on export, Brender will generate an obj file for each object in the scene snapshotted at every time step. These can be used as desired, but we usually import them into blender for nicer looking renders than what tour basi OpenGL settup provides.

**NOTE**: The export directory must already exist in the file system, the simulation will not generate a directory for you.

## Contact
If you would like to contact us for anything regarding EoL-Cloth free to email us.
If you have any code specific comments or find any bugs, please specifically contact Nick Weidner via [GitHub](https://github.com/weidnern "Nick Weidner GitHub") or Email

