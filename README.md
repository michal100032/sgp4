## sgp4
A simple implementation of the SGP4 perturbation model in modern C++ as a programming exercise. The goal is to make the original C-style code by David Vallando more readable but not less accurate.

The project still is under development, although time conversion utilites and TLE parsing are already working*

\* kind of
  
If you have any suggestions or questions feel free to post an issue or email me.

### Build
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

### Dependencies
Currently none
