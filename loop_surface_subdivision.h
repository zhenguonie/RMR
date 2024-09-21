#pragma once
#pragma once
//#include <strategy.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <map>

class Loop
{
private:
	typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;
	MyMesh mesh;
	std::map<int, OpenMesh::Vec3d> facePoints;
	std::map<int, OpenMesh::Vec3d> edgePoints;
	std::map<int, OpenMesh::Vec3d> vertexPoints;
	int times;
public:
	Loop(Data* data_)
};