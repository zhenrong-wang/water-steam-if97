# Water and Steam Thermophysical Property Calculation

# 1. Background

The thermophysical properties of water and steam are very important to researchers and engineers in various domains (energy, power, chemistry, materials ...).

Usually, researchers and engineers rely on the [NIST Refprop](https://www.nist.gov/srd/refprop) to calculate the properties. However, Refprop is commercial and not open-source. Therefore, I wrote this code to provide an alternative to the Refprop.

# 2. Brief Intro

**Program Name**: Ultra-high-resolution Water and Steam Thermophysical Property Calculation

**Purpose**: Calculate the **very-high-resolution** thermophysical properties at given point(s)

**License**: MIT

**Technical Reference**: [IAPWS IF97](http://www.iapws.org/)

# 3. How-To: Build, Run, and Use

## 3.1 Build

### 3.1.1 Prerequisites

You need a C compiler to build. 

- For Microsoft Windows users, [mingw](https://sourceforge.net/projects/mingw/) is a good choice
- For GNU/Linux Distro or other *nix users, the [GNU Compiler Collections](https://gcc.gnu.org/), known as gcc, is a perfect one
- For macOS users, [clang](https://clang.llvm.org/) is easy to install and use (brew is not needed to install clang on macOS).

### 3.1.2 Build Guide

1. Use `git` to clone this code: `git clone https://github.com/zhenrong-wang/water-steam-if97.git`
2. Build command example: `gcc *.c -o my-if97.exe -lm`

Note: the `-lm` may not be valid for Windows or macOS. It is necessary for GNU/Linux distros.

## 3.2 Run

### An Example for UNIX-like OS:

Suppose the working directory is `/home/amber/water-steam/`

Suppose the absolute path of the executable is `/home/amber/bin/my-if97.exe`

- `$ cd /home/amber/water-steam/`
- `$ vim _input.dat` # Create an input file named '**_input.dat**' with some input data (See Section 3.3 below)
- `$ /home/amber/bin/my-if97.exe`

### An Example for Windows:

Suppose the working directory is `C:\Users\amber\water-steam\`

Suppose the absolute path of the executable is `C:\Users\amber\bin\my-if97.exe`

- Open a Command Prompt Window
- `cd C:\Users\amber\water-steam\`
- `notepad _input.dat`
- `C:\Users\amber\bin\my-if97.exe`

## 3.3 Use

For the **_input.dat** file, you need to input the data in lines. A single line refers to a point. Strict format:

`POINT_TYPE(int),PARAM1(float/double),PARAM2 (float/double)` 

**Please separate the params by a comma ','**

POINT_TYPE: 

- 1 Pressure and Temperature
- 2 Pressure and Density
- 3 Pressure and Specific Internal Energy
- 4 Pressure and Specific Enthalpy
- 5 Pressure and Specific Entropy
- 6 Temperature and Density
- 7 Temperature and Specific Internal Energy
- 8 Temperature and Specific Enthalpy
- 9 Temperature and Specific Entropy
- 10 Specific Enthalpy and Specific Entropy
- 11 Pressure and Steam Dryness
- 12 Temperature and Steam Dryness

Example: `1,100000,300`

Example: `1,1e5,3e2`

Units:

- Pressure(p): Pa
- Temperature(t): K
- Density(r): kg/m3
- Specific Internal Energy(u): J/kg
- Specific Enthalpy(h): J/kg
- Specific Entropy(s): J/(kg.K)

## 3.4 Output

The program writes out a file named '**_properties.dat**' in the working directory. 

# 4 Bugs and Communications

This code was written years ago, and I didn't continuously develop it. Therefore, bugs may occur. 

If you are interested in this project, please submit issues to this repo. I'd be glad to communicate on any issues.

Or, you can also email me: zhenrongwang@live.com

