# Subdivision Shading for Catmull-Clark and Loop Subdivision Surfaces with Semi-Sharp Creases

https://www.mdpi.com/2073-431X/12/4/85

## How to use this project:
1. Use Qt to run the program.
2. Load a model using the “Load OBJ file” button.
3. In case the OBJ file didn’t specify sharpness values, the “Load sharpness file” button can be used to load a sharpness preset in the form of a .sharp file (see specification below).
4. Use the “Geometric subdivision” field to subdivide model on screen.
4. Enable normal subdivision with the options in the “Subdivision left view” field.
5. Specify the preferred blending method with the options in the “Blend normals” field.
6. Visualize certain properties of your subdivided model with the options under “Draw mode” and “Colour mode”.

## Specification of .sharp files:
The sharp files are a preset of sharpness values for a model. The file consists of a list of sharpness values for each face in the model in the order the faces appear in their respective .OBJ file. Each line consists of either of the following 2 options:
 - “fs x” where x is a set of digits specifying the sharpness values of a face in clockwise order. For example: “fs 1 2 3 2 5”.
- “as y x” where x is the same as before and y is a positive integer. This will assign the sharpness values x to the next y faces.
If not all faces are covered by the sharp file the remaining sharpness values in the model will be set to 0. If 2 faces specify different values for a vertex they share, the highest value will be used.

## Specification of .OBJ file:
OBJ files use their standard layout, with one addition: A sharpness preset can be appended at the bottom of the file using the same rules as used in a .sharp file
