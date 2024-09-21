#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO/OFF.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Subdivision_method_3/subdivision_methods_3.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO/PLY.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Surface_mesh_deformation.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <boost/iterator/function_output_iterator.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/copyleft/cgal/RemeshSelfIntersectionsParam.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/opengl/glfw/Viewer.h>

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <ctime>
#include <omp.h>

unsigned int flag_saved = 0;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3     Point;
typedef K::Vector_3    Vector;
typedef K::Aff_transformation_3 Transformation;

typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::Search_traits_3<K> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits> Neighbor_search;
typedef Neighbor_search::Tree NeiTree;
typedef CGAL::Kd_tree<Traits> KdTree;


typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor            face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor        halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor            edge_descriptor;
typedef Mesh::Vertex_index Vertex_index;

namespace PMP = CGAL::Polygon_mesh_processing;

struct halfedge2edge
{
    halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges)
        : m_mesh(m), m_edges(edges)
    {}
    void operator()(const halfedge_descriptor& h) const
    {
        m_edges.push_back(edge(h, m_mesh));
    }
    const Mesh& m_mesh;
    std::vector<edge_descriptor>& m_edges;
};

void saveData(std::string fileName, Eigen::MatrixXd matrix)
{
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",", "\n");

    std::ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}

void SAVE_P_AND_N(Mesh& mesh)
{
    unsigned int N = mesh.number_of_vertices();
    auto vnormals = mesh.add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
    PMP::compute_vertex_normals(mesh, vnormals);
    unsigned int op = 0;
    Eigen::MatrixX3d Out_points(N, 3);
    Eigen::MatrixX3d pt_normal(N, 3);
    for (auto v : mesh.vertices()) {
        Out_points.row(op) << mesh.point(v).x(), mesh.point(v).y(), mesh.point(v).z();
        pt_normal.row(op) << vnormals[v].x(), vnormals[v].y(), vnormals[v].z();
        op++;
    }
    saveData("Out_points.csv", Out_points);
    saveData("mesh_normal.csv", pt_normal);
}

struct Coff_Matrix
{
    unsigned int n;
    SpMat LK;
    Eigen::MatrixX3d B;
    Eigen::VectorXd prop;
    Eigen::MatrixXd pt_normal;
};

Mesh Init_Evolution(Vector site, const double R, const unsigned int iter = 3)
{
    //Using a regular icosahedron to approximate a sphere
    Mesh mesh;
    const double m = sqrt(50 - 10 * sqrt(5)) / 10 * R;
    const double n = sqrt(50 + 10 * sqrt(5)) / 10 * R;
    std::vector<Point> points = {
        Point(m,0,n),
        Point(m,0,-n),
        Point(-m,0,n),
        Point(-m,0,-n),
        Point(0,n,m),
        Point(0,-n,m),
        Point(0,n,-m),
        Point(0,-n,-m),
        Point(n,m,0),
        Point(-n,m,0),
        Point(n,-m,0),
        Point(-n,-m,0)
    };
    std::vector<std::vector<std::size_t>> faces = {
        {6, 4, 8}, {9, 4, 6}, {6, 3, 9}, {6, 1, 3}, {6, 8, 1},
        {8, 10, 1}, {8, 0, 10}, {8, 4, 0}, {4, 2, 0}, {4, 9, 2},
        {9, 11, 2}, {9, 3, 11}, {3, 1, 7}, {1, 10, 7}, {10, 0, 5},
        {0, 2, 5}, {2, 11, 5}, {3, 7, 11}, {5, 11, 7}, {10, 5, 7}
    };
    for (const auto& p : points)
    {
        auto& pm = p;
        mesh.add_vertex(pm + site);
    }
    for (const auto& f : faces)
        mesh.add_face(Mesh::Vertex_index(f[0]), Mesh::Vertex_index(f[1]), Mesh::Vertex_index(f[2]));
    for (int i = 0; i < iter; ++i) {
        CGAL::Subdivision_method_3::Loop_subdivision(mesh);
    }
    return mesh;
}

void Igl_Show(Mesh& mesh)
{

    Eigen::MatrixXd V, V_p;
    Eigen::MatrixXi F, F_p;
    const std::string file_saved = "Cube.off";
    if (!igl::readPLY("PointClouds\\Box.ply", V_p, F_p))
    {
        std::cerr << "Failed to open the file!" << std::endl;
        return;
    }
    if (!CGAL::write_off(mesh, file_saved))
    {
        std::cout << "Write OFF file Error!" << std::endl;
        return;
    }
    igl::readOFF(file_saved, V, F);
    //// Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_points(V_p, Eigen::RowVector3d(1, 0, 0));
    viewer.data().point_size = 1.0;
    viewer.launch();
}

void HolesFilling(Mesh& mesh)
{
    unsigned int nb_holes = 0;
    for (halfedge_descriptor h : halfedges(mesh))
    {
        if (CGAL::is_border(h, mesh))
        {
            std::vector<face_descriptor>  patch_facets;
            std::vector<vertex_descriptor> patch_vertices;
            bool success = std::get<0>(PMP::triangulate_refine_and_fair_hole(mesh, h,
                std::back_inserter(patch_facets),
                std::back_inserter(patch_vertices),
                CGAL::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)).geom_traits(K())));

            CGAL_assertion(CGAL::is_valid_polygon_mesh(mesh));

            std::cout << "* FILL HOLE NUMBER " << ++nb_holes << std::endl;
            std::cout << "  Number of facets in constructed patch: " << patch_facets.size() << std::endl;
            std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
            std::cout << "  Is fairing successful: " << success << std::endl;
        }
    }

    CGAL_assertion(CGAL::is_valid_polygon_mesh(mesh));
    std::cout << std::endl;
    std::cout << nb_holes << " holes have been filled" << std::endl;

}

void Border_Stitching(Mesh& mesh)
{
    std::vector<halfedge_descriptor> border_cycles;
    std::vector<halfedge_descriptor> border_i;
    std::vector<std::vector<halfedge_descriptor>> boundary;
    //Collect one halfedge per boundary cycle
    PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
    int boundary_count = 0;
    for (halfedge_descriptor h : border_cycles) {
        halfedge_descriptor cur = h;
        do {
            border_i.push_back(cur);
            std::cout << "Halfedge " << cur << std::endl;
            cur = mesh.next(cur);
        } while (cur != h);
        boundary.push_back(border_i);
        border_i.clear();
    }
    std::vector<halfedge_descriptor> b1;
    std::vector<halfedge_descriptor> b2;
    //Judge whether the number of vertices on two boundaries is equal
    if (boundary[0].size() <= boundary[1].size())
    {
        //std::cout << "Size is same!!!" << std::endl;
        b1 = boundary[0];
        b2 = boundary[1];

    }
    else
    {
        //std::cout << "Size is not same!!!" << std::endl;
        b1 = boundary[1];
        b2 = boundary[0];
    }
    std::vector<halfedge_descriptor> b2_sorted;
    //Find the nearest boundary pair
    Point s1 = mesh.point(mesh.source(b1[0]));
    halfedge_descriptor closet_halfedge;
    double min_distance = std::numeric_limits<double>::max();
    for (halfedge_descriptor h2 : b2)
    {
        Point t2 = mesh.point(mesh.target(h2));
        double dist = CGAL::squared_distance(s1, t2);
        if (dist < min_distance)
        {
            min_distance = dist;
            closet_halfedge = h2;
        }
    }
    halfedge_descriptor cur_ = closet_halfedge;
    do {
        b2_sorted.push_back(cur_);
        cur_ = mesh.prev(cur_);
    } while (cur_ != closet_halfedge);
    //Stiching
    for (unsigned int hi = 0; hi < b1.size(); hi++)
    {
        //1 - Obtain the two vertices (starting and ending points) of the half chain
        Vertex_index src1 = mesh.source(b1[hi]);
        Vertex_index tgt1 = mesh.target(b1[hi]);

        Vertex_index src2 = mesh.source(b2_sorted[hi]);
        Vertex_index tgt2 = mesh.target(b2_sorted[hi]);

        //2 - Generate two triangular facets between these four vertices
        mesh.add_face(src1, tgt1, src2);
        mesh.add_face(src1, src2, tgt2);

    }
}

double Mixed_Distance(Neighbor_search& search, const unsigned int K, const double lamda = 0.75)
{
    //Calculate the mixed distance field
    double result = 0.0, Wpq_sum = 0.0;
    std::vector<double> vec;
    for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    {
        vec.push_back(std::sqrt(it->second));
    }
    for (std::vector<double>::iterator it = vec.begin(); it != vec.end(); it++)
    {
        Wpq_sum += 1 / pow(*it, 2);
    }

    for (std::vector<double>::iterator it = vec.begin(); it != vec.end(); it++)
    {
        result += ((1 / pow(*it, 2)) / Wpq_sum) * lamda * (*it);
    }
    return result;
}

double Activate_Distance(double dis, double down_dis = 0)
{
    if (dis < down_dis)
    {
        return  0;
    }
    else
    {
        return dis * 0.9;
    }
}

void Isotropic_Remeshing(Mesh& mesh, const double target_edge_length = 1, const unsigned int nb_iter = 5)
{
    std::cout << "Split border...";
    std::vector<edge_descriptor> border;
    PMP::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(halfedge2edge(mesh, border)));
    PMP::split_long_edges(border, target_edge_length, mesh);
    std::cout << "done." << std::endl;
    std::cout << "Start remeshing : " << " (" << num_faces(mesh) << " faces)..." << std::endl;
    PMP::isotropic_remeshing(faces(mesh), target_edge_length, mesh,
        CGAL::parameters::number_of_iterations(nb_iter)
        .protect_constraints(true)); //i.e. protect border, here
    std::cout << "Remeshing done." << std::endl;
    return;
}

void Mesh_Smoothing(Mesh& mesh, const unsigned int nb_iterations = 1, const double time = 0.3)
{
    std::set<Mesh::Vertex_index> constrained_vertices;
    for (Mesh::Vertex_index v : vertices(mesh))
    {
        if (is_border(v, mesh))
            constrained_vertices.insert(v);
    }
    std::cout << "Constraining: " << constrained_vertices.size() << " border vertices" << std::endl;
    CGAL::Boolean_property_map<std::set<Mesh::Vertex_index> > vcmap(constrained_vertices);
    std::cout << "Smoothing shape... (" << nb_iterations << " iterations)" << std::endl;
    PMP::smooth_shape(mesh, time, CGAL::parameters::number_of_iterations(nb_iterations)
        .vertex_is_constrained_map(vcmap));
    CGAL::IO::write_polygon_mesh("mesh_shape_smoothed.off", mesh, CGAL::parameters::stream_precision(17));
    std::cout << "Done!" << std::endl;
}

double compute_cotangent(const Vector& u, const Vector& v)
{
    // Compute dot product of u and v
    double dot_product = u * v;
    // Compute cross product of u and v
    Vector cross_product = CGAL::cross_product(u, v);
    // Compute the norm of the cross product
    double cross_product_norm = std::sqrt(cross_product.squared_length());
    // Compute cotangent value
    double cotangent = dot_product / cross_product_norm;
    return cotangent;
}

Coff_Matrix Phase1(Mesh& mesh, const Point_set& PCD, double wf, double ws, double wp, double eps = 1, double w_mixdistance = 1)
{
    //Parameter initialization and variable space allocation
    unsigned int N = mesh.number_of_vertices();//Number of evolution source vertices
    unsigned int nx = N, ny = 3 * N;//Matrix pre allocated space size
    unsigned int K_NUM = 10;// K Neighbor points
    unsigned int iik = 0;//Parallel computing parameters
    //
    SpMat LK(nx, ny);
    std::vector<T> tripletList;
    Eigen::MatrixX3d Ind(N, 3);
    Eigen::MatrixX3d Bf(N, 3), Bs(N, 3), Bp(N, 3);
    Eigen::MatrixXd pt_normal(N, 3);
    //
    tripletList.reserve(ny);
    LK.reserve(ny);
    Ind.setZero();
    Bf.setZero(); Bs.setZero(); Bp.setZero();
    auto vnormals = mesh.add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;

    //KD - Tree : Store input point cloud
    NeiTree tree(PCD.points().begin(), PCD.points().end());
    PMP::compute_vertex_normals(mesh, vnormals);

    //Parallel computing lock
    std::vector<vertex_descriptor> VD;
    for (auto v : mesh.vertices())
    {
        VD.push_back(v);
    }
#pragma omp parallel num_threads(4)
    {
#pragma omp for reduction(+:iik)    
        // Evolutionary source vertex traversal
        for (int ii = 0; ii < N; ii++)
        {
            // Locate vertex information
            auto v = VD[ii];
            Point point_before = mesh.point(v);
            Vector pvi(point_before.x(), point_before.y(), point_before.z());

            // Bp
            Bp(ii, 0) = point_before[0]; Bp(ii, 1) = point_before[1]; Bp(ii, 2) = point_before[2];

            // Adaptive mixed distance field
            Eigen::Vector3d P_search, nn;
            nn << vnormals[v].x(), vnormals[v].y(), vnormals[v].z();
            Neighbor_search search(tree, point_before, K_NUM);
            int signal = 1;
            for (Neighbor_search::iterator it = search.begin(); it != search.end(); it++) {
                //std::cout << " Point " << it->first << std::endl;
                P_search << it->first.x(), it->first.y(), it->first.z();
                if (P_search.dot(nn) < 0)
                {
                    signal = -1;
                    break;
                }
            }
            double mixed_dis = Mixed_Distance(search, K_NUM, w_mixdistance);
            if (mixed_dis < eps)
            {
                if (signal < 0)
                {
                    mixed_dis = signal * 0.8 * mixed_dis;
                }
                else {
                    mixed_dis = 0;
                }
            }

            // Bf
            Vector Evolution_step(pvi + (mixed_dis)*vnormals[v]);
            pt_normal.row(ii) << vnormals[v].x(), vnormals[v].y(), vnormals[v].z();
            Bf(ii, 0) = Evolution_step[0]; Bf(ii, 1) = Evolution_step[1]; Bf(ii, 2) = Evolution_step[2];

            // Bs
            std::vector<vertex_descriptor> one_ring;
            std::vector<double> w_cot;
            double Asi = 0;
            Vector Bsi(0, 0, 0);
            for (auto vj : mesh.vertices_around_target(mesh.halfedge(v)))
            {
                one_ring.push_back(vj);
            }
            unsigned int num_1_ring = one_ring.size();
            one_ring.insert(one_ring.begin(), one_ring.back());
            one_ring.push_back(one_ring[1]);
            for (unsigned int i = 1; i < one_ring.size() - 1; i++)
            {
                Vector ua = point_before - mesh.point(one_ring[i - 1]);
                Vector va = mesh.point(one_ring[i]) - mesh.point(one_ring[i - 1]);
                double cota = compute_cotangent(ua, va);
                Vector ub = point_before - mesh.point(one_ring[i + 1]);
                Vector vb = mesh.point(one_ring[i]) - mesh.point(one_ring[i + 1]);
                double cotb = compute_cotangent(ub, vb);
                Asi += (cota + cotb) / 2;
                Bsi += (cota + cotb) * (mesh.point(one_ring[i]) - Point(0, 0, 0)) / 2;
            }
            unsigned int Dii = one_ring.size() - 2;

            // Coefficient assignment
            double r = pow(wf, 2) + pow(ws * Asi / Dii, 2) + pow(wp, 2);
            Ind(ii, 0) = wf / r;
            Ind(ii, 1) = ws * Asi / (Dii * r);
            Ind(ii, 2) = wp / r;
            Bs(ii, 0) = Bsi[0] / Dii; Bs(ii, 1) = Bsi[1] / Dii; Bs(ii, 2) = Bsi[2] / Dii;

            one_ring.clear();
            iik++;
        }
    }

    //Parallel computing determines whether data is lost
    if (iik != N)
    {
        std::cerr << "iik=" << iik << "  and  N=" << N << std::endl;
        throw std::runtime_error("------Error: iik and N not equal-------");
    }

    // Coefficient sorting
    mesh.remove_property_map(vnormals);
    for (unsigned int iit = 0; iit < N; iit++)
    {
        tripletList.push_back(T(iit, iit, Ind(iit, 0)));
        tripletList.push_back(T(iit, iit + N, Ind(iit, 1)));
        tripletList.push_back(T(iit, iit + 2 * N, Ind(iit, 2)));
    }
    LK.setFromTriplets(tripletList.begin(), tripletList.end());
    Eigen::MatrixX3d B(3 * N, 3);
    B << wf * Bf, ws* Bs, wp* Bp;

    Coff_Matrix result;
    result.n = N;
    result.pt_normal = pt_normal;

    // LK * V = B  where V = LK-1 * B is the closed-form solution
    result.LK = LK;
    result.B = B;

    return result;
}

Eigen::MatrixX3d Phase2(Mesh& mesh, const Point_set& PCD, double r = 5, double sq_control = 500.0)
{
    //Parameter initialization and variable space allocation
    unsigned int N = mesh.number_of_vertices();
    unsigned int iik = 0;
    Eigen::MatrixX3d Out_points(N, 3);
    Eigen::MatrixX3d pt_xyz(PCD.points().size(), 3);
    Eigen::MatrixX3d pt_normal(N, 3);
    Out_points.setZero();
    auto vnormals = mesh.add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
    PMP::compute_vertex_normals(mesh, vnormals);
    unsigned int ci = 0;
    for (auto it = PCD.points().begin(); it != PCD.points().end(); ++it)
    {
        Point point = *it;
        pt_xyz.row(ci) << point.x(), point.y(), point.z();
        ++ci;
    }
    std::vector<vertex_descriptor> VD;
    for (auto v : mesh.vertices())
    {
        VD.push_back(v);
    }
#pragma omp parallel num_threads(4)
    {
#pragma omp for reduction(+:iik)
        for (int ii = 0; ii < N; ii++)
        {
            // Locate vertex information
            auto v = VD[ii];
            Point point_before = mesh.point(v);
            Vector pvi(point_before.x(), point_before.y(), point_before.z());

            double aver_dis = 0;
            Eigen::VectorXd P0(3);
            Eigen::Vector3d n;
            P0 << point_before.x(), point_before.y(), point_before.z();
            n << vnormals[v].x(), vnormals[v].y(), vnormals[v].z();
            pt_normal.row(ii) << vnormals[v].x(), vnormals[v].y(), vnormals[v].z();

            //Retain point clouds with acute angles to the normal direction of the point cloud
            Eigen::MatrixXd PQ = pt_xyz.rowwise() - P0.transpose();
            Eigen::VectorXd direction = PQ * n;
            int num_p_p = ((direction.array() > 0).cast<int>()).sum();
            Eigen::MatrixXd New_PQ(num_p_p, 3);
            unsigned int i_newpq = 0;
            for (unsigned int i_posi = 0; i_posi < direction.size(); i_posi++)
            {
                if (direction(i_posi) > 0)
                {
                    New_PQ.row(i_newpq) << PQ(i_posi, 0), PQ(i_posi, 1), PQ(i_posi, 2);
                    i_newpq++;
                }
            }

            //Calculate the set of point clouds with a distance less than r from the normal line
            Eigen::VectorXd distance(num_p_p);
            for (unsigned int i_row = 0; i_row < num_p_p; i_row++)
            {
                Eigen::Vector3d temp_pq = New_PQ.row(i_row);
                distance(i_row) = temp_pq.cross(n).norm();
            }
            i_newpq = 0;
            int num_p_d = ((distance.array() < r).cast<int>()).sum();
            if (num_p_d == 0)
            {
                aver_dis = 0;
            }
            else
            {
                Eigen::MatrixXd New_ptxyz(num_p_d, 3);
                Eigen::VectorXd distance0(num_p_d);
                distance0.setZero();
                for (unsigned int i_dis = 0; i_dis < distance.size(); i_dis++)
                {
                    if (distance(i_dis) < r)
                    {
                        New_ptxyz.row(i_newpq) << New_PQ(i_dis, 0), New_PQ(i_dis, 1), New_PQ(i_dis, 2);
                        distance0(i_newpq) = New_ptxyz.row(i_newpq).norm();
                        i_newpq++;
                    }
                    else {

                    }
                }
                //Controlling evolutionary distance through variance
                Eigen::VectorXd diff_dis = distance0.array() - distance0.mean();
                double sq_dis = diff_dis.array().square().sum() / num_p_d;
                if (sq_dis > sq_control) {
                    aver_dis = distance0.minCoeff();
                }
                else {
                    aver_dis = distance0.sum() / num_p_d;
                    aver_dis = Activate_Distance(aver_dis);
                }
            }
            //Projection Evolution
            Vector p_after(pvi + (aver_dis)*vnormals[v]);
            Out_points(ii, 0) = p_after[0];
            Out_points(ii, 1) = p_after[1];
            Out_points(ii, 2) = p_after[2];
            iik++;
        }
    }
    return Out_points;
}


int main()
{
    // Initialize evolutionary source
    Vector site = Vector(0, 0, 0);  //Center coordinates
    int R = 10.0;  //Sphere radius
    Mesh mesh = Init_Evolution(site, R);

    // Import input point cloud
    Point_set PCD;
    std::ifstream f(CGAL::data_file_path("PointClouds\\Box.ply"), std::ios_base::binary);

    if (!CGAL::IO::read_PLY(f, PCD))
    {
        std::cerr << "Can't read input file ply" << std::endl;
        return EXIT_FAILURE;
    }

    // Initialize parameters
    // 1- The number of iterations in two stages
    unsigned int ITER1 = 20, ITER2 = 1;
    unsigned int STOP_ITER = ITER1 + ITER2;
    unsigned int RECORD_ITER = STOP_ITER;
    // 2- Cut off accuracy
    double eps = 5;
    double eps_all = 0.3;
    // 3- Regulatory factor
    double wf = 0.3, ws = 0.9, wp = 0.8;
    // 4- The length of mesh edge
    double l_control = 5;
    // Iteration
    while (STOP_ITER--)
    {
        // Info.
        std::cout << "Cu - Iteraton - " << STOP_ITER << "- Left - " << std::endl;

        unsigned int N = mesh.number_of_vertices();
        Eigen::MatrixX3d V(N, 3);
        Eigen::MatrixX3d B(3 * N, 3);
        if (STOP_ITER >= ITER2)
        {
            Coff_Matrix AB;
            try {
                AB = Phase1(mesh, PCD, wf, ws, wp, eps);
            }
            catch (const std::runtime_error& e) {
                std::cerr << "Caught an exception: " << e.what() << std::endl;
            }

            SpMat LK = AB.LK;   B = AB.B;
            if (N == AB.n)
            {
                V = LK * B;
                unsigned int i_v = 0;
                double l_all = 0;
                for (auto v : mesh.vertices()) {
                    Point point_new(V(i_v, 0), V(i_v, 1), V(i_v, 2));
                    if (1)
                    {
                        Vector delta_d = point_new - mesh.point(v);
                        l_all += delta_d.squared_length();
                    }
                    mesh.point(v) = point_new;
                    i_v++;
                }
                std::cout << "Error = " << l_all / i_v << std::endl;
            }
        }
        else
        {
            unsigned int i_v = 0;
            V = Phase2(mesh, PCD);
            int tes = mesh.vertices().size();
            for (auto v : mesh.vertices()) {
                Point point_new(V(i_v, 0), V(i_v, 1), V(i_v, 2));
                mesh.point(v) = point_new;
                i_v++;
            }
        }
        //Topology Refining
        if (STOP_ITER % 2 == 0 || STOP_ITER <= ITER2) {
            bool intersecting = PMP::does_self_intersect(mesh);
            if (intersecting)
            {
                std::cout << "There are self-intersections." << std::endl;
                std::vector<std::pair<face_descriptor, face_descriptor>> intersected_tris;
                std::set<face_descriptor> faces_to_remove;
                PMP::self_intersections(mesh, std::back_inserter(intersected_tris));
                for (const auto& p : intersected_tris)
                {
                    auto hf1 = mesh.halfedge(p.first);
                    auto hf2 = mesh.halfedge(p.second);
                    auto fc_start1 = mesh.faces_around_face(hf1);
                    auto fc_start2 = mesh.faces_around_face(hf2);
                    for (auto fc1 : fc_start1)
                    {
                        faces_to_remove.emplace(fc1);
                    }
                    for (auto fc2 : fc_start2)
                    {
                        faces_to_remove.emplace(fc2);
                    }
                }
                for (auto& f : faces_to_remove) {
                    CGAL::Euler::remove_face(mesh.halfedge(f), mesh);
                }

                Mesh::Property_map<face_descriptor, std::size_t>fccmap =
                    mesh.add_property_map<face_descriptor, std::size_t>("f:CC").first;
                std::size_t num = PMP::connected_components(mesh, fccmap);

                std::cout << "number of connected components:" << num << std::endl;
                std::cerr << "- We keep the largest components" << std::endl;
                PMP::keep_largest_connected_components(mesh, 1);
                if (STOP_ITER != 0)
                {
                    //Border_Stitching(mesh);
                }
                HolesFilling(mesh);
                Isotropic_Remeshing(mesh, l_control, 3);
            }
            else
            {
                std::cout << "There are NO self-intersections." << std::endl;
                Isotropic_Remeshing(mesh, l_control, 3);
            }
        }

    }
    //Surface topology optimization
    Isotropic_Remeshing(mesh, l_control / 1.5, 3);
    Mesh_Smoothing(mesh);

    //SAVE_P_AND_N(mesh);
    Igl_Show(mesh);
    return 0;
}
