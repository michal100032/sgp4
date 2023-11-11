## sgp4
A simple implementation of the SGP4 perturbation model in modern C++ as a programming exercise. The goal is to make the original C-style code by David Vallando more readable but not less accurate.

The project still is under development, although time conversion utilites and TLE parsing are already working*

\* kind of
  
If you have any suggestions or questions feel free to post an issue or email me.

### Clone and build
First clone the repository
```
git clone git@github.com:michal100032/sgp4.git
cd sgp4
```
Then, run cmake
```
cmake -B out/build
```
And finally build the project with make
```
cd out/build
make
```
To run the test program just type
```
tests/sgp4_test
```
Alternatively to build the project you can simply open it in Visual Studio.

### Sources
Sources that I base my code on
* A series of tutorials on orbital propagation on Celestrak (http://celestrak.org/columns/)
* Spacetrak Report no. 3 (https://celestrak.org/NORAD/documentation/spacetrk.pdf). This is where the SGP4 along with SDP4 (the deep space version) were described. Unfortunatelly the algorithm for generating the TLEs has been modified over the years and thus the original algorithms are not entirely valid anymore. The report contains FORTRAN implementation of the algorithms.
* Revisiting Spacetrack Report no. 3 (http://celestrak.org/publications/AIAA/2006-6753/). This is an update of the algorithms. There is also a C++ implementation last updated in 202 
### Dependencies
Currently none
