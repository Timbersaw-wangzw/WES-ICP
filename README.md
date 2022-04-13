# WES-ICP
This repository includes the source code of the paper Robust Point Clouds Registration with Point-to-point $l_p$ Distance Constraints in Large-scale Metrology.
This code was written by MATLAB.
# Introduction and Usage
The program includes the following three source point clouds and the same target point clouds.
Those source point clouds have been coarse transformed, but they have different overlapping areas between target point clouds.
```
source point clouds:
- coarse_source_points0.70.txt
- coarse_source_points0.6.txt
- coarse_source_points0.5.txt
target point clouds:
- target points.txt
```
The directory `github_repo` is the lie algebra library.
The program includes three algorithms
```
1. sparse point to point
2. sparse point to plane
3. ours(WES-ICP)
```

Our **WES-ICP** can restrain the sliding and escape the local minima. Specifically, **sparse point to point** is easily trapped into local minima and **sparse point-to-plane** slides along the tangent space. Those phenomenons can be seen in our project.
you can run directly `mainICP.m` to see those phenomenons.
# Phenomenons
![inital position](https://github.com/Timbersaw-wangzw/WES-ICP-M/blob/master/inital%20pose.jpg)
![slding](https://github.com/Timbersaw-wangzw/WES-ICP-M/blob/master/sparse%20point%20to%20plane%20gap.jpg)

![local minima](https://github.com/Timbersaw-wangzw/WES-ICP-M/blob/master/sparse%20point%20to%20point.jpg)

![restrain the sliding and escape the local minima](https://github.com/Timbersaw-wangzw/WES-ICP-M/blob/master/WES-ICP.jpg)
