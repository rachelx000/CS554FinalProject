# CS 554 Final Project

## Adaptive Local Remeshing for Smooth Silhouette on 3D Triangle Meshes

### Rachel Xing

<table>
   <td width="50%">
      <img src="figure/bunny_smooth_sil.gif" alt="bunny_smooth_sil" style="width: 100%; float: left; margin-right: 10px;" />
   </td>
   <td width="50%">
      <img src="figure/feline_smooth_sil.gif" alt="feline_smooth_sil" style="width: 100%; float: left; margin-right: 10px;" />
   </td>
</table>

&nbsp;
### How to compile the program
This is based on the graphics program called **learnply** for scientific visualization provided by Prof. Eugene Zhang @ Oregon State University.

CMake is used to compile the program, so please make sure a version of 3.10 or higher is installed.

The general steps to build the "learnply" program are as follows:
1. Download this project as a `.zip` file and unzip it
2. Direct to the folder for the unzipped file in the terminal
3. In the terminal, run the following command to compile the program and then run it with Visual Studio on Windows or in the terminal on Mac/Linux:
#### Windows (with Visual Studio)
```bat
md build
cd build
cmake -A Win32 ..
```
**Note:** You must include the `-A Win32` switch in the CMake command, as the executable must be a 32-bit executable.

Then, open the compiled `learnply.sln` using Visual Studio. When building and running on Visual Studio, be sure to select the **learnply** project 
as the **startup project** through the context menu by right clicking on it in the solution explorer.

#### Mac/Linux
```bash
mkdir build
cd build
cmake ..
make
```
To run the compiled program, run the following command in the terminal:
```bash
./learnply
```

### How to run the program with different models
Open `learnply.cpp` in the folder `learnply` using any code editor, and go to the line 105:
```C++
this_file = fopen("../tempmodels/xxx.ply", "r");
```
Replace the `xxx` with the name of a model in the folder `tempmodels`. Then, compile and run the project as described above.
The folder includes both coarse (`coarse_xxx.ply`) and fine (`fine_xxx.ply`) versions of torus, Stanford bunny, happy Budda, 
feline, and dragon I created using MeshLab. You can also add any valid `.ply` model to `tempmodels` and use it in the same way.

### How to use the program

1. **Silhouette Smoothing Mode (this project)**: press the key `m` to toggle the adaptive local remeshing pipeline for silhouette smoothing. Use the following keys:
   1. key `1` - Render shaded mesh after remeshing
   2. key `2` - Render wireframe after remeshing (newly added geometry highlighted in red)
   3. key `3` - Combined shaded and wireframe rendering
   4. key `s` - Toggle rendering of silhouette lines generated by the pipeline. You can hide the model by pressing any other number key and show the drawn silhouette only.
   5. key `x` - Toggle anti-aliasing for current rendering
2. **Default Rendering mode**: Visualization of the input model without silhouette smoothing. Use the following keys:
   1. key `1` - Shaded rendering of the original mesh
   2. key `2` - Wireframe rendering of the original mesh
   3. key `3` - Polygon ID map visualization
   4. key `4` - Barycentric coordinates map visualization
   5. key `5` - Normal map (face normals) visualization
   6. key `6` - Apply 3D Checkerboard texture using L-sized cubes
      - Default `L = 1.0`
      - Press `p` to double L (increase cube size)
      - Press `l` to halve L (decrease cube size)
   7. key `7` - Visualize Gaussian curvature computed by Angle Deficit method
   8. key `8` - Visualize Gaussian curvature computed by the Valence Deficit method
   8. key `s` - Toggle rendering of silhouette lines found by face-based silhouette extraction algorithm. You can hide the model by pressing any other number key and show the drawn silhouette only.
   9. key `x` - Toggle anti-aliasing for current rendering



