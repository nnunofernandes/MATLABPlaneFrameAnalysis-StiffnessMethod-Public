# Plane Frame Analysis Using MATLAB

Teaching materials for solving two-dimensional frames with the **direct stiffness method**. The repository combines worked exercises, MATLAB starter and solution scripts, guided tutorials, and a PDF explaining the numerical method alongside independent analytical solutions.

Both examples include **axial deformation and bending**, with three degrees of freedom at each node: horizontal displacement, vertical displacement, and rotation.

## Contents

- [Learning objectives](#learning-objectives)
- [Teaching materials](#teaching-materials)
- [Getting started](#getting-started)
- [Exercise 16-2](#exercise-16-2)
- [Exercise 16-16](#exercise-16-16)
- [Stiffness method workflow](#stiffness-method-workflow)
- [Units and sign conventions](#units-and-sign-conventions)
- [Outputs and verification](#outputs-and-verification)
- [Classroom activities](#classroom-activities)
- [Common mistakes](#common-mistakes)

## Learning objectives

After completing the exercises, students should be able to:

1. Define node coordinates, element connectivity, and global DOF numbering.
2. Form the local stiffness matrix of an Euler–Bernoulli plane-frame element.
3. Transform stiffness matrices and loads into global coordinates.
4. Convert distributed loads into consistent nodal forces **and moments**.
5. Assemble the global equations and apply support conditions.
6. Solve nodal motions and recover reactions and member-end actions.
7. Interpret deformed shapes and check physical equilibrium.
8. Compare numerical results with analytical calculations.

## Teaching materials

| Material | Purpose |
| --- | --- |
| Tutorial PDF | Step-by-step explanations, equations, diagrams, and analytical solutions for both exercises. |
| Guided MATLAB tutorial | Sectioned code explaining each computational stage and its physical meaning. |
| MATLAB before solution | Starter script for students to complete. |
| MATLAB after solution | Completed implementation for checking the exercise. |

### Suggested file organization

The paths below define the intended repository organization. Place the teaching files under these names when publishing, or update these paths to match the final filenames.

| Path | Material |
| --- | --- |
| `README.md` | Repository overview and exercise statements. |
| `docs/MATLAB_Plane_Frame_Step_by_Step_Guide.pdf` | Shared tutorial PDF. |
| `tutorials/exercise_16_2_plane_frame.m` | Guided tutorial for Exercise 16-2. |
| `tutorials/exercise_16_16_plane_frame.m` | Guided tutorial for Exercise 16-16. |
| `exercises/exercise_16_2_before.m` | Exercise 16-2 starter script. |
| `exercises/exercise_16_2_after.m` | Exercise 16-2 completed solution. |
| `exercises/exercise_16_16_before.m` | Exercise 16-16 starter script. |
| `exercises/exercise_16_16_after.m` | Exercise 16-16 completed solution. |

The guided tutorials use `%%` sections and explanatory comments in ordinary `.m` files. A native MATLAB Live Script edition, if included, uses the `.mlx` extension.

## Getting started

Use MATLAB with support for local functions in scripts. The supplied guided `.m` tutorials include their own supporting functions and require no additional toolboxes.

1. Download or clone the repository.
2. Open MATLAB and set the repository root as the current folder.
3. Read the relevant exercise statement below.
4. Attempt its **before** script.
5. Work through the guided tutorial, running the `%%` sections **in order**.
6. Compare your work with the **after** script and the analytical explanation in the PDF.

For the suggested folder layout, run either guided tutorial from the repository root:

```matlab
run(fullfile('tutorials', 'exercise_16_2_plane_frame.m'))
```

```matlab
run(fullfile('tutorials', 'exercise_16_16_plane_frame.m'))
```

Run one exercise at a time. Each guided script initializes its workspace and figure state.

### Display options

```matlab
showMatrices = true;   % Display intermediate matrices and load vectors
makePlots = true;      % Plot the model and deformed shape
deformationScale = 0;  % Choose deformation magnification automatically
```

Set `showMatrices = false` for shorter output. Magnification affects only the plot; displacement tables report the actual calculated motions.

## Exercise 16-2

### Displacement of a cantilever frame

Determine the **horizontal displacement, vertical displacement, and rotation at node 1** of an L-shaped frame.

The horizontal beam extends **3 m** from node 1 to node 2. A **4 m** column connects node 2 to the fixed base at node 3. Node 2 is a rigid corner. At node 1, apply an **18 kN force to the left** and a **27 kN force downward**.

All members have:

| Property | Value |
| --- | --- |
| Young’s modulus, $E$ | $200\ \mathrm{GPa}$ |
| Cross-sectional area, $A$ | $12.9\times10^3\ \mathrm{mm^2}$ |
| Second moment of area, $I$ | $300\times10^6\ \mathrm{mm^4}$ |

### Model definition

Coordinates are in millimetres. The DOF lists follow `[ux uy theta]`.

| Node | x | y | Global DOFs | Condition |
| --- | ---: | ---: | --- | --- |
| 1 | 0 | 4000 | `[1 2 3]` | Loaded free end |
| 2 | 3000 | 4000 | `[4 5 6]` | Rigid corner |
| 3 | 3000 | 0 | `[7 8 9]` | Fixed support |

| Element | Connectivity | Length | Local x direction |
| --- | --- | ---: | --- |
| 1 | Node 1 → Node 2 | 3000 mm | Right |
| 2 | Node 2 → Node 3 | 4000 mm | Down |

**Restrained DOFs:** `[7 8 9]`.  
**Free DOFs:** `[1 2 3 4 5 6]`.

### Tasks

1. Define the geometry, material properties, connectivity, and DOF map.
2. Form and transform both element stiffness matrices.
3. Assemble the nine-DOF system and apply the fixed-base restraints.
4. Solve the six unknown motions and report the three components at node 1.
5. Recover the base reactions and local member-end actions.
6. Verify equilibrium and compare with the analytical solution.

### Analytical parallel

Let $P=18000\ \mathrm{N}$ and $Q=27000\ \mathrm{N}$ be the positive magnitudes of the leftward and downward loads, with $B=3000\ \mathrm{mm}$ and $H=4000\ \mathrm{mm}$.

The corner rotation and tip motions are:

$$
\theta_2=\frac{PH^2}{2EI}+\frac{QBH}{EI},
$$

$$
u_{x1}=-\frac{PH^3}{3EI}-\frac{QBH^2}{2EI}-\frac{PB}{EA},
$$

$$
u_{y1}=-\frac{QH}{EA}-B\theta_2-\frac{QB^3}{3EI},
\qquad
\theta_1=\theta_2+\frac{QB^2}{2EI}.
$$

These terms describe column bending, axial shortening, corner rotation, and beam bending. The PDF also derives the translations from strain energy using Castigliano’s theorem.

<details>
<summary><strong>Reference results for Exercise 16-2</strong></summary>

Positive directions are right, up, and counterclockwise.

| Requested quantity | Result | Direction |
| --- | ---: | --- |
| $u_{x1}$ | −17.22093023 mm | Left |
| $u_{y1}$ | −27.49186047 mm | Down |
| $\theta_1$ | +0.009825000 rad | Counterclockwise |

| Support | Horizontal reaction | Vertical reaction | Moment |
| --- | ---: | ---: | ---: |
| Node 3 | +18.000000 kN | +27.000000 kN | −153.000000 kN m |

The corner rotation contributes $-B\theta_2=-23.4\ \mathrm{mm}$ to the vertical tip motion.

</details>

## Exercise 16-16

### Reactions of a portal frame

Determine the **horizontal and vertical reactions at supports 1 and 4**.

The portal has two **3 m columns** joined by a **2 m beam**. The bases at nodes 1 and 4 are **pins**. The upper connections at nodes 2 and 3 are **rigid**. A uniform horizontal load of **6 kN/m to the right** acts along the entire left column.

All members have:

| Property | Value |
| --- | --- |
| Young’s modulus, $E$ | $200\ \mathrm{GPa}$ |
| Cross-sectional area, $A$ | $9.68\times10^3\ \mathrm{mm^2}$ |
| Second moment of area, $I$ | $300\times10^6\ \mathrm{mm^4}$ |

### Model definition

Coordinates are in millimetres. The nonsequential DOF numbering follows the supplied exercise figure.

| Node | x | y | Global DOFs `[ux uy theta]` | Condition |
| --- | ---: | ---: | --- | --- |
| 1 | 0 | 0 | `[11 12 8]` | Left pin support |
| 2 | 0 | 3000 | `[1 2 3]` | Upper left rigid joint |
| 3 | 2000 | 3000 | `[4 5 6]` | Upper right rigid joint |
| 4 | 2000 | 0 | `[9 10 7]` | Right pin support |

| Element | Connectivity | Length | Local x direction |
| --- | --- | ---: | --- |
| 1 | Node 1 → Node 2 | 3000 mm | Up |
| 2 | Node 2 → Node 3 | 2000 mm | Right |
| 3 | Node 3 → Node 4 | 3000 mm | Down |

**Restrained DOFs:** `[9 10 11 12]`.  
**Free DOFs:** `[1 2 3 4 5 6 7 8]`.

Base rotations **7 and 8 remain free**. Prescribing them as zero would model fixed bases.

### Distributed loading

Element 1 points upward, so its local positive y direction points left. The rightward load therefore requires a negative local transverse component:

```matlab
w = 6;                     % N/mm = 6 kN/m
qLocal = [0 -w;             % Left column
          0  0;
          0  0];
```

The consistent nodal load vector gives horizontal forces of **+9 kN at each end**, with an equivalent moment of **−4.5 kN m at the lower end** and **+4.5 kN m at the upper end**.

These are equivalent applied loads. The recovered **physical moments at both pin supports are zero**.

### Tasks

1. Define the model while preserving the exercise DOF numbering.
2. Restrain base translations and retain the unknown base rotations.
3. Convert the member load into consistent nodal forces and moments.
4. Assemble the twelve-DOF system and solve the eight unknown motions.
5. Recover reactions using the complete global load vector.
6. Recover local member-end actions and check the pin moments.
7. Verify whole-frame equilibrium and compare with the force method.

### Analytical parallel

Release the horizontal restraint at node 4 and choose $X=R_{x4}$ as the redundant, positive right. With $H=3000\ \mathrm{mm}$, $B=2000\ \mathrm{mm}$, and $w=6\ \mathrm{N/mm}$, the unit-load method gives:

$$
\delta_0=\frac{w}{EI}\left(\frac{5H^4}{24}+\frac{H^3B}{4}\right),
\qquad
f_{xx}=\frac{2H^3/3+H^2B}{EI}+\frac{B}{EA}.
$$

Restoring the support requires zero horizontal displacement:

$$
\delta_0+f_{xx}X=0,
\qquad
X=-\frac{\delta_0}{f_{xx}}.
$$

The PDF develops the bending-moment products, flexibility integrals, and remaining reactions. The term $B/EA$ retains beam axial flexibility.

<details>
<summary><strong>Reference results for Exercise 16-16</strong></summary>

Positive directions are right, up, and counterclockwise.

| Support | Horizontal reaction | Vertical reaction | Moment |
| --- | ---: | ---: | ---: |
| Node 1 | −12.94620144 kN | −13.50000000 kN | 0 |
| Node 4 | −5.05379856 kN | +13.50000000 kN | 0 |

Both horizontal reactions act left. The negative vertical reaction at node 1 means the support must hold the frame down.

The compatibility calculation gives $\delta_0=3.037500\ \mathrm{mm}$ and $f_{xx}=0.000601033058\ \mathrm{mm/N}$.

</details>

## Stiffness method workflow

The guided tutorials follow the same computational sequence:

1. Define material and cross-section data.
2. Enter node coordinates and member connectivity.
3. Map nodal motions to global DOFs.
4. Specify support motions and external loads.
5. Calculate element lengths, orientations, and DOF lists.
6. Form local stiffness matrices and transform them to global axes.
7. Construct consistent nodal vectors for member loads.
8. Assemble the global stiffness matrix and load vector.
9. Partition the equations and solve the free displacements.
10. Recover reactions and local member-end actions.
11. Check equilibrium and compare with the analytical solution.
12. Plot and interpret the deformed shape.

For $\mathbf{d}_{\ell}=\mathbf{T}\mathbf{d}_g$:

$$
\mathbf{k}_g=\mathbf{T}^{T}\mathbf{k}_{\ell}\mathbf{T},
\qquad
\mathbf{f}_g=\mathbf{T}^{T}\mathbf{f}_{\ell}.
$$

The reduced system, reactions, and member-end actions are:

$$
\mathbf{K}_{ff}\mathbf{d}_f=\mathbf{F}_f-\mathbf{K}_{fc}\mathbf{d}_c,
$$

$$
\mathbf{R}=\mathbf{K}\mathbf{d}-\mathbf{F},
\qquad
\mathbf{q}_{\mathrm{end}}=\mathbf{k}_{\ell}\mathbf{d}_{\ell}-\mathbf{f}_{\ell}.
$$

MATLAB solves the reduced system with backslash:

```matlab
d(freeDofs) = Kff \ rhs;
```

### Modeling assumptions

- Straight, prismatic Euler–Bernoulli members in a two-dimensional frame.
- Linear elastic material behavior and small displacements and rotations.
- Axial and bending deformation included.
- Rigid member-to-member connections at the upper corners.
- Shear deformation, geometric nonlinearity, and member self-weight excluded.

## Units and sign conventions

| Quantity | Internal unit or convention |
| --- | --- |
| Length and translation | mm |
| Force | N |
| Distributed force | N/mm |
| Moment | N mm |
| Rotation | rad |
| Young’s modulus | N/mm² |
| Positive global translation | Right and up |
| Positive rotation and moment | Counterclockwise |
| Local x direction | From element node i to node j |
| Local y direction | 90° counterclockwise from local x |
| Local DOF order | `[ui vi thetai uj vj thetaj]` |
| Local end-action order | `[Ni Vi Mi Nj Vj Mj]` |

The scripts retain the exercise **DOF numbers** and explicitly use counterclockwise-positive rotations. Local end actions act **on the member**; their signs differ from some conventions for internal force diagrams.

## Outputs and verification

| Output | Description |
| --- | --- |
| `displacementTable` | Nodal translations in mm and rotations in rad. |
| `reactionTable` | Support forces in kN and moments in kN m. |
| `endForceTable` | Local axial forces, shears, and moments at both member ends. |
| `K`, `F`, `Kff` | Full and reduced stiffness-system data. |
| `d`, `R` | Complete displacement vector and reaction/residual vector. |
| `elementData` | Element geometry, transformations, stiffness matrices, and loads. |
| Frame plots | Original geometry and magnified curved deformation. |

Verification includes matrix symmetry, positive definiteness of the restrained stiffness, equilibrium at free DOFs, whole-frame force and moment balance, and an analytical comparison. The portal also checks that both pin-end moments are zero.

The formulations have been independently checked using numerical integration of element stiffness and loads, cantilever and energy calculations, and the force method. Direct MATLAB runtime verification remains to be performed.

Small floating-point residuals near zero are expected. Printed digits help compare implementations; they do not imply corresponding precision in real material or geometric data.

## Classroom activities

1. Sketch the likely deformation before running either script.
2. Explain why the vertical load in Exercise 16-2 also creates horizontal motion.
3. Double $E$ and explain why displacements halve while these reaction values remain unchanged.
4. Increase $A$ while keeping $I$ fixed and inspect which results change.
5. Explain why the portal load needs equivalent nodal moments.
6. Explain why a pin’s rotation remains unknown when its physical moment is zero.
7. Identify the shared DOFs receiving stiffness contributions from multiple members.
8. Reproduce an analytical displacement or reaction without using the completed script.

When modifying loads or geometry, update the analytical-check inputs consistently before comparing results.

## Common mistakes

| Mistake | Correction |
| --- | --- |
| Mixing metres and millimetres | Convert all inputs before assembly. |
| Constraining portal base rotations | Keep rotational DOFs 7 and 8 free. |
| Incorrect local load direction | The upward left column requires `qy = -6` for a rightward load. |
| Omitting equivalent nodal moments | Include both end-moment terms for the uniform transverse load. |
| Removing loads at restrained DOFs | Keep them in the full `F` used to recover `R`. |
| Overwriting shared-joint stiffness | Add each member matrix into `K(edofs,edofs)`. |
| Using `kLocal*dLocal` alone for a loaded member | Subtract `fLocal` to recover physical end actions. |
| Reading DOF order as node order | Gather results using `dofMap`. |
| Measuring actual deflection from the magnified plot | Use the numerical displacement table. |

---

Educational material for structural mechanics, matrix structural analysis, and introductory finite element methods.
