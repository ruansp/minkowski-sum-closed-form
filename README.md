# Closed-Form Minkowski Sums of Convex Bodies with Smooth Positively Curved Boundaries
[![View Closed-form Minkowski sums on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/116000-closed-form-minkowski-sums)

Matlab implementation for exact closed-form Minkowski sums of general convex bodies with smooth positively curved boundaries. Scripts are provided for all the results shown in the article.

__The article__: [Publisher](https://www.sciencedirect.com/science/article/abs/pii/S0010448521001445), [Arxiv](https://arxiv.org/abs/2012.15461)

## Authors:
+ Sipu Ruan (repository maintainer): National University of Singapore, ruansp@nus.edu.sg
+ Greogory S. Chirikjian: National University of Singapore, mpegre@nus.edu.sg

## Introduction
This work derives closed-form parametric formulas for the Minkowski sums of convex bodies in $d$-dimensional Euclidean space with boundaries that are smooth and have all positive sectional curvatures at every point. Under these conditions, there is a unique relationship between the position of each boundary point and the surface normal. 

The main results are presented as two theorems:
+ Theorem 1 (Theorem 4.1 in the article): directly parameterizes the Minkowski sums using the unit normal vector at each surface point. Although simple to express mathematically, such a parameterization is not always practical to obtain computationally;
+ Theorem 2 (Theorem 4.3 in the article): derives a more useful parametric closed-form expression using the gradient that is not normalized. 

Demonstrations of the proposed closed-form Minkowski sums:
<table>
  <tr>
    <td><img src="/misc/demo_minksum_thm1.png" alt="Demonstration of normal-parameterized closed-form Minkowski sums" width="300"/></td>
    <td><img src="/misc/demo_minksum_thm2.png" alt="Demonstration of gradient-parameterized closed-form Minkowski sums" width="300"/></td>
  </tr>
  <tr>
    <th>Normal parameterization</th>
    <th>Gradient parameterization</th>
  </tr>
</table>


In the special case of two ellipsoids, the proposed expressions are identical to those derived previously using geometric interpretations. In order to examine the results, numerical validations and comparisons of the Minkowski sums between two superquadric bodies are conducted. 

Two applications are discussed and demonstrated:
+ Generate configuration space obstacles in motion planning problems;
+ Improve performance for optimization-based collision detection algorithms.

Demonstrations of collision detection using the proposed closed-form Minkowski sums:
<table>
  <tr>
    <td><img src="/misc/app_collision_3D_demo_mink_ray.png" alt="Collision detetion using closed-form Minkowski sums and ray-casting method" width="300"/></td>
    <td><img src="/misc/app_collision_3D_demo_mink_normal.png" alt="Collision detetion using closed-form Minkowski sums and common normal concept" width="300"/></td>
  </tr>
  <tr>
    <th>Ray-casting method</th>
    <th>Common-normal concept</th>
  </tr>
</table>

## Implementations
### Classes and functions (located in "/src" folder)
+ __Ellipse.m__, __Ellipsoid.m__, __EllipsoidND.m__: Class that defines ellipsoids in 2D, 3D and N-D
+ __SuperEllipse.m__, __SuperQuadrics.m__: Class that defines superquadrics in 2D and 3D
+ __MinkSumClosedForm.m__: Class that defines the closed-form Minkowski sums between two convex bodies with smooth and positively curved boundaries (works for any dimension)
+ __MinkSumDefinition.m__: Function that computes Minkowski sums using convex hull method (works for any dimension)
+ __MinkSumEdgeSort2D.m__: Function that computes Minkowski sums using edge sorting method (works only in 2D)

### Scripts (located in "/test" folder)
#### Demonstrations 
+ __demo_minksum_definition.m__: Demonstrations of different algorithms to compute Minkowski sums (convex hull, edge sort and closed-form) using 2D superellipses *[Produce raw plots for Figure 1]*
+ __demo_minksum_theorem.m__: Demonstrations of Theorem 4.1 and 4.3 using 2D superellipses *[Produce raw plots for Figure 2]*
+ __demo_minksum_sq_2D.m__, __demo_minksum_sq_3D.m__: Demonstrations for the closed-form Minkowski sums between two superellipses in 2D/3D *[Produce Figure 3 for 2D case and Figure 4 (a)-(b) for 3D case]*
+ __demo_minksum_primitives_3D.m__: Demonstrations for the closed-form Minkowski sums between two basic geometric primitives approximated by superquadrics in 3D *[Produce Figure 4 (c)-(h)]*

#### Verifications and Benchmarks
+ __test_geometry_fitting_error.m__: Comparisons of fitting errors using different geometric models, i.e., superquadrics, convex polyhedra and geometric primitives *[Reproduce results in Section 6.2]*
+ __verify_minksum_2D.m__, __verify_minksum_3D.m__: Verifications of closed-form Minkowski sums in 2D/3D based on averaged distance from implicit surface (DI) and contact number (Nc), i.e., Equation (36) *[Reproduce results in Section 6.1, Eq. (36)]*
+ __verify_minksum_2D_kissing_point.m__, __verify_minksum_3D_kissing_point.m__: Verifications of closed-form Minkowski sums in 2D/3D based on the error of kissing point, i.e., Equation (37) *[Reproduce results in Section 6.1, Eq. (37)]*
+ __bench_mink2D.m__, __bench_mink3D.m__: Benchmark between different algorithms of Minkowski sums for 2D/3D cases *[Produce Figures 5-7 for ellipsoid-elliosoid (EE), ellipsoid-superquadric (ES) and superquadric-superquadric (SQ) cases]*

#### Applications
+ __app_c_obstacle_generate_SE2.m__, __app_c_obstacle_generate_SE3.m__: Application of C-obstacle generation in 2D/3D cases *[Produce Figure 8, reproduce results in Table 2]*
+ __app_collision_3D_demo.m__: Demonstration of collision detection between two superquadrics using Minkowski sum-based methods and common normal concept *[Produce Figure 9]*
+ __app_collision_3D.m__: Benchmarks of different collision detection algorithms between two superquadrics, i.e., Minkowski sum-based, implicit surface, common normal, GJK, algebraic separation condition (specific for two ellipsoids) *[Reproduce results in Table 3]*

## Citation
If you find our work interesting or have used in your research, please consider citing our article:

- Ruan, S. and Chirikjian, G.S., 2021. Closed-Form Minkowski Sums of Convex Bodies with Smooth Positively Curved Boundaries. Computer-Aided Design, p.103133.

- BibTex:
```
@article{ruan2021closed,
  title={Closed-Form Minkowski Sums of Convex Bodies with Smooth Positively Curved Boundaries},
  author={Ruan, Sipu and Chirikjian, Gregory S},
  journal={Computer-Aided Design},
  pages={103133},
  year={2021},
  publisher={Elsevier}
}
```
