# A Robust Surface Mesh Reconstruction Method for Real-Time 3D Visualization in Minimally Invasive Surgery: Addressing Noise, Topology, and Deformation Accuracy

## Abstract
Lightweight triangular mesh model has great potential for real-time 3D visualization of the lesion during minimally invasive surgery(MIS). However, the blurred tissue boundaries, high imaging noise, and unoriented points of medical images seriously affect the accuracy and topological quality of surface reconstruction, which will lead to inaccurate lesion localisation. In this paper, we develop a robust and high-topology-quality triangular mesh reconstruction method that provides a good deformable expression model for intraoperative real-time 3D visualization. This method first approximates the model prototype under the guidance of an unsigned distance field by simulating inflation. And then, we design a variance controlled cylindrical domain projection search (VC-CDPS) method to complete the final surface fitting. In addition, we incorporate a topology optimization in the iterative reconstruction process to ensure smoothness and good topology of the reconstruction model. We validated the method on the geometric model with high noise and the human organ model manually segmented by novice doctors. The results show that the reconstructed model has better surface quality and anti-robustness. Moreover, we design a comparison experiment of model deformation and propose a metric to measure the topological quality of the model. Through tissue experiments in vitro, we explore the relationship between the topological quality and the accuracy of deformation. The results show that the deformation accuracy is positively correlated with the topological quality.

![](https://github.com/Scalpelapex/Images/blob/main/VC_CDPS/Overview.jpg)

## Environment setup

Clone the repo: 
> https://github.com/zhenguonie/VFCM-CHPS

Open this project with Visual Studio 2019 on Windows x64

## Libraries

Eigen >= 3.4.0

Libigl >= 2.4.0

CGAL >= 5.5.2

## Datasets
The model dataset provided in this project is composed of multiple parts:
1) The partial-liver model and the kidney model are extracted from the dataset AMOS:
> http://ncov-ai.big.ac.cn/download

2) The partial-liver model and the kidney model are extracted from the dataset CC-CCII:
> https://zenodo.org/records/7155725#.Y0OOCOxBztM

3)The sphere, cube and  toroidal cube are created by us. The sphere model and cube model have high noise along the direction of thickness. The toroidal cube model is sparse and have a genu.

**The point clouds extraction version of the dataset is saved in folder PointClouds**

## Comments
We have provided five test cases, and their corresponding code and parameters are saved in the following five CPP files：

Case01_Liver.cpp

Case02_Lung.cpp

Case03_Cube.cpp

Case04_Sphere.cpp

Case05_Square_Torus.cpp


The input and output of the model can be modified in the code.

## Results

### Process Diagram

![](https://github.com/Scalpelapex/Images/blob/main/VC_CDPS/Surface.gif)

### Comparison Results

![](https://github.com/Scalpelapex/Images/blob/main/VC_CDPS/Results.jpg)

## Citation

If you found this code helpful, please consider citing:

```

```
