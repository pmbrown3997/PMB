1.) Main Script: MHD_model.m
This script utilizes a combination of finite differences to calculate the velocities from the momentum equation and numerical integration to calculate the pressure from the force balance equation.

2.) Advection.m: calculates the necessary advection operations involved in the continuity equation in terms of velocities in the x, y directions and densities of the plasma.

3.) HHD.m: Helmholtz Hodges Decomposition decouples the pressure and velocity fields. If the computation preserves the orthogonality of the decomposition, any error in one of the terms is not reflected in the other.