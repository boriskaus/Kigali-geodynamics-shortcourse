# Building more complicated LaMEM geometries
*Boris Kaus*

## 1. Introduction
Until now, we focussed on relatively simple model geometries where we can use geometrical objects such as layers, spheres, hexahedrons etc. Yet, in many cases it would be nice to have more flexibility in creating model setups. What LaMEM needs is information at every gridpoint what the rock `Phase` and `Temperature` is. 
This can be done with the `GeophysicalModelGenerator.jl` julia package.


## 2. Subduction setup commented

The full julia file to create the julia setup for the simulation is the file `CreateMarkers_Subduction_Linear_FreeSlip_parallel.jl`. 

Initially, we load the package and indicate the LaMEM parameter file:

```julia
using GeophysicalModelGenerator  

# Load LaMEM particles grid
ParamFile   =   "Subduction2D_FreeSlip_Particles_Linear_DirectSolver.dat"
Grid        =   ReadLaMEM_InputFile(ParamFile)
```
`Grid` is a structure thatr contains information about the LaMEM grid, based on the parameters in the LaMEM ParamFile. 

Next, we specify some useful parameters for this setup:
```julia
# Specify slab parameters
Trench_x_location   = -500;     # trench location
Length_Subduct_Slab =  200;     # length of subducted slab
Length_Horiz_Slab   =  1500;    # length of overriding plate of slab
Width_Slab          =  750;     # Width of slab (in case we run a 3D model)         

SubductionAngle     =   34;     # Subduction angle
ThicknessCrust      =   10;     
ThicknessML         =   75;     # Thickness of mantle lithosphere
T_mantle            =   1350;   # in Celcius
T_surface           =   0;

ThicknessSlab = ThicknessCrust+ThicknessML
```

Now we initialize the matrixes that contain the phases and temperature info
```julia
Phases      =   zeros(Int64, size(Grid.X));             # Rock numbers
Temp        =   ones(Float64,size(Grid.X))*T_mantle;    # Temperature in C    
```

In this setup, we create a slab by using the `AddBox!` function which is part of GMG. See the GMG help page on how to use it. It makes it a bit easier to generate a lithosphere with layers of different thickness
```julia
# Create horizontal part of slab with crust & mantle lithosphere
AddBox!(Phases,Temp,Grid,
        xlim=(Trench_x_location, Trench_x_location+Length_Horiz_Slab), 
        zlim=(-ThicknessSlab   , 0.0),
        phase=LithosphericPhases(Layers=[ThicknessCrust ThicknessML], Phases=[1 2 0]) );               

# Add inclined part of slab                            
AddBox!(Phases,Temp,Grid,
        xlim=(Trench_x_location-Length_Subduct_Slab, Trench_x_location), 
        zlim=(-ThicknessSlab   , 0.0),
        DipAngle=-SubductionAngle,
        Origin=(Trench_x_location,0,0),
        phase=LithosphericPhases(Layers=[ThicknessCrust ThicknessML], Phases=[1 2 0]) );               
```

In order to create a 3D model setup, we need to create a 3D Cartesian Data structure that contains `Phases` and `Temp` as fields:
```julia
# Save julia setup 
Model3D     =   CartData(Grid, (Phases=Phases,Temp=Temp))   # Create LaMEM model:
Write_Paraview(Model3D,"LaMEM_ModelSetup")                  # Save model to paraview   (load with opening LaMEM_ModelSetup.vts in paraview)  
```

Finally, we save the LaMEM markers to disk. On one processor, that is easy and can be done with one line:
```julia
# Save LaMEM markers
Save_LaMEMMarkersParallel(Model3D)                          # Create LaMEM marker input on 1 core
```


If you, however, want to run the simulation in parallel, the situation is a bit more complicated as every processor needs to read its own portion of the particles. You do this by creating a partitioning file `PartFile` and passing that along. 
```julia
julia> PartFile = run_lamem_save_grid(ParamFile,4)
julia> Save_LaMEMMarkersParallel(Model3D, PartitioningFile=PartFile)

Writing LaMEM marker file -> ./markers/mdb.00000000.dat
Writing LaMEM marker file -> ./markers/mdb.00000001.dat
Writing LaMEM marker file -> ./markers/mdb.00000002.dat
Writing LaMEM marker file -> ./markers/mdb.00000003.dat
```

### 2.1 Subduction exercise
Run the above setup on 4 or 8 cores and change the subduction angle.
