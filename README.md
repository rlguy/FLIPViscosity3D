# FLIPViscosity3D

This program is a basic free-surface liquid FLIP fluid simulator with viscosity. The fluid simulation program exports the simulation particle positions as a sequence of .OBJ or .PLY meshes containing only vertices. The particles can then be meshed and rendered in a separate program.

This program implements the following methods for the variational pressure solver and the variational viscosity solver:

* [Fast Variational Framework for Accurate Solid-Fluid Coupling](https://hal.archives-ouvertes.fr/file/index/docid/384725/filename/variationalFluids.pdf) by Batty, Bertails, and Bridson
* [Accurate Viscous Free Surfaces for Buckling, Coiling, and Rotating Liquids](https://cs.uwaterloo.ca/~c2batty/papers/BattyBridson08.pdf) by Batty and Bridson

Much of the simulation method was adapted from code available from [Christopher Batty's GitHub page](https://github.com/christopherbatty) and also code available from my [GridFluidSim3D](https://github.com/rlguy/GridFluidSim3D) project.

## Gallery

The following images were generated from particle data output by the program and rendered using [Blender](http://www.blender.org).

[![alt tag](http://i.imgur.com/EOgftEL.jpg)](http://i.imgur.com/R97LLa5.jpg) [![alt tag](http://i.imgur.com/K5OLqDj.jpg)](http://i.imgur.com/weVF42W.jpg)

[![alt tag](http://i.imgur.com/vivBTG0.jpg)](http://i.imgur.com/M2G2tm4.jpg) [![alt tag](http://i.imgur.com/7DnmSfm.jpg)](http://i.imgur.com/lkaoatu.jpg)

An animation of a buckling honey simulation can be viewed [here](https://www.youtube.com/watch?v=Oxsr4m-s3C8).

## Features
Below is a list of features implemented in the simulator.

* FLIP (Fluid-Implicit-Particle) simulation method
* Accurate variational pressure solve for free-surfaces and curved boundaries
* Accurate viscosity for buckling, coiling, and rotating liquids
* Initialize fluids and solid boundaries from triangle meshes

## Dependencies

Everything required to run this program is included in the repository. A compiler that supports C++11 is required to build the program.

## Installation

This program uses the [CMake](https://cmake.org/) utility to generate the appropriate solution, project, or Makefiles for your system. A Makefile is also provided.

Sample triangle meshes are located in the ```sample_meshes``` directory.

The default simulation will drop a mass of fluid in the shape of the Stanford Bunny inside of a spherical container.

[![alt tag](http://i.imgur.com/sRZi5bQ.jpg)](http://i.imgur.com/rRRnJXs.jpg)

## Rendering in Blender

An example script for how to import the particle meshes into [Blender](http://blender.org) for rendering is located [here](src/blender/render_particles.py). This script will import a .ply or .obj mesh into Blender and duplicate a sphere over the particles.

Usage: load the script into the Blender text editor, edit the ```MESH_DIRECTORY``` variable to point to the directory containing the simulation meshes, and press the 'Run Script' button.