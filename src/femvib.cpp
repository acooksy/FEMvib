// femvib.cpp is the main FEMvib program, constructing the A and B matrices from the PES and G-matrix 
//  elements, then diagonalizing with the SLEPc eigensolver.
//
// Original author: Dong Xu, 2007
// Extensive modifications: Peter Zajac, 2013
// Rewritten with updated calls for compatibility with current libraries: ALC, 2024
// 13-Aug-2025 v2.2.00 ALC
// Added support for target energy, rather than always calculating eigenvalues from ground state on up.
//  Adds a new input option "-w targ_en" where targ_en = target energy of the SLEPc eigensolver.
//  In general this means nev eigenpairs will be obtained, roughly centered on targ_en.
// Since the above revision changed the workflow to leave libMesh earlier and let PETSc do all the
//  post-processing, the output format has changed.  In particular, eigenstates are now printed
//  directly to ascii xyz files, rather than as ExodusII files that required an additional translator script.
//
// The original FEMvib C++ code femvib.cpp was built on top of:
// <h1>Eigenproblems Example 2 - Solving a generalized Eigen Problem</h1>
// \author Steffen Petersen
// \date 2006
// In this example some eigenvalues for a generalized symmetric
// eigenvalue problem A*x=lambda*B*x are computed, where the
// matrices A and B are assembled according to stiffness and
// mass matrix, respectively.
//
//  FEMvib should be run using the femvib_run.sh script from a directory with a subdirectory "./results"
//  which houses these files:
//    input  (values for each of the input options)
//    coordinates  (the geometries at all points on the PES)
//    energyP  (the energies at all points on the PES)
//  Formats for these last two files are given below or see examples in tests/:
//
//  coordinates:
//    (# geometries)  (# atoms in molecule)  (# PES coordinates = PES dimensionality
//    [1st geometry is a reference geometry (typically a stationary point) in xyz format with atomic masses:]
//       mass_1  x_1  y_1  z_1
//       mass_2  x_2  y_2  z_2
//         ...
//       mass_N  x_N  y_N  z_N
//    (value coord 1)  (value coord 2)  (value coord 3)      [assuming a 3D PES]
//    [2nd geometry is a geometry from some point on PES in xyz format with no leading column]
//       x_1  y_1  z_1
//       x_2  y_2  z_2
//         ...
//       x_N  y_N  z_N
//    [Repeat for all geometries.  
//     The reference geometry should be given again, so it can be linked to the correct PES coordinates.
//     The total number of lines in the coordinates file is therefore (N_atom + 1) * (N_points + 1)
//     where N_atom is the number of atoms in the molecule and N_points is the number of points in the PES.]
//
//  energyP:
//      (# energies)
//      (value coord 1)  (value coord 2)  (value coord 3)  (relative energy 1 in cm-1)
//      (value coord 1)  (value coord 2)  (value coord 3)  (relative energy 2 in cm-1)
//         ...
//      (value coord 1)  (value coord 2)  (value coord 3)  (relative energy N in cm-1)
//    [where the relative energy is typically computed relative to the equilibrium geometry.]
//     The total number of lines in the energyP file is therefore (N_points + 1) ]
//

// libMesh include files.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/eigen_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/slepc_eigen_solver.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/enum_eigen_solver_type.h"
#include "libmesh/getpot.h"
#include "libmesh/point.h"
#include "libmesh/point_locator_base.h"
// the above all from the original eigenproblems_ex2 file
//// the below added for femvib
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <getopt.h>
#include <chrono>
#include <ctime>
//#include "Eigen/src/SparseCore/SparseView.h"   // this one breaks with an ISO C++ error
#include "Eigen/src/Householder/HouseholderSequence.h"
#include "Eigen/src/Core/VectorwiseOp.h"
#include "Eigen/src/Core/util/Memory.h"
#include "Eigen/src/Core/util/Macros.h"
#include "Eigen/src/Core/CoreEvaluators.h"
#include "include/nr3_alc.h"
#include "include/ludcmp.h"
//#include "svd.h"
//#include "helpf.h"
#include "include/pointSelector2_alc.h"
//#include "imls.h"
#include "include/interp_rbf.h"
#include "include/krig.h"
//#include "moment.h"
//#include "dataman.h"
#include "include/duplicates_alc.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/slepc_macro.h"
#include "libmesh/elem.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

std::unordered_map<std::string, std::string> variables;

// Function prototype.  This is the function that will assemble
// the eigen system. Here, we will simply assemble a mass matrix.
void assemble_mass(EquationSystems & es,
                   const std::string & system_name);

void write_petsc_eigenvector_ascii(const EPS& eps,
                                   PetscInt eig_index,
                                   const libMesh::EquationSystems& equation_systems,
                                   const std::string& system_name,
                                   const std::string& var_name,
                                   const std::string& output_filename);

void parseFile(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    
    while (std::getline(file, line)) {
        size_t pos = line.find('=');
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos);
            std::string value = line.substr(pos + 1);
            variables[key] = value;
        }
    }
}


int main (int argc, char ** argv)
{
  int ndim;
  auto start = std::chrono::system_clock::now();  // start timekeeping

  // Initialize libMesh and the dependent libraries.
  LibMeshInit init (argc, argv);

#ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
  // SLEPc currently gives us a nasty crash with Real==float
  libmesh_example_requires(false, "--disable-singleprecision");
#endif

  // Initiate log file
std::ofstream logfile("results/femvib.log", std::ios::app);
std::string versioninfo = "FEMvib v2.2.00\n 13-Aug-2025 \n";
std::string authorlist = "Authors: Dong Xu, Peter Zajac, Jernej Stare, Andrew L. Cooksy\n";
std::string citation = "Citation: Comp. Phys. Commun. 180, 2079-2094 (2009). (doi:10.1016/j.cpc.2009.06.010.)\n";
std::string contactinfo = "Contact: acooksy@sdsu.edu\n";
std::string licenseinfo = "Released under the GNU GENERAL PUBLIC LICENSE version 3 2007\n\n";
logfile << versioninfo << authorlist << citation << contactinfo << licenseinfo
              << std::endl;
std::time_t start_time = std::chrono::system_clock::to_time_t(start);
logfile << "Starting computation at " << std::ctime(&start_time) << "\n";
std::cout << "Starting computation at " << std::ctime(&start_time) << "\n";

  // Tell the user what we are doing.
  libMesh::out << "Running " << argv[0];

  for (int i=1; i<argc; i++)
    libMesh::out << " " << argv[i];

  libMesh::out << std::endl << std::endl;

  // READ ARGUMENTS
  // Read:
  // -ndim [# dimensions]
  // -n [# eigenvalues to be computed]
  // -nx -ny -nz  [# points in mesh along each coord]
  // -targ_en  [target energy eigenvalue]
  // -x0 -x1 -y0 -y1 -z0 -z1  [lower / upper bounds along each coord]
  // -pbcx -pbcy -pbcz  [flags for periodfic boundary conditions along each coord]
  // -elemTypeStr  [FEM element type]
  // -npr  [# wavefunction files to print; default is n]
  // -npr1  [1st wavefunction to print]; output will be for npr
  //             eigenstates starting from npr1

  int opt;
  std::string inputFile;
  ndim = 2;
  std::string elemTypeStr = "QUAD9";
  int nev = 5;
  int npr = nev;
  int npr1 = 0;
  int nx = 10;
  int ny = 10;
  int nz = 10;
  double targ_en = 0.0;
  double x0 = 0.0;
  double x1 = 1.0;
  double y0 = 0.0;
  double y1 = 1.0;
  double z0 = 0.0;
  double z1 = 1.0;
  int pbc_x = 0;
  int pbc_y = 0;
  int pbc_z = 0;

// switch statements should use single-character flags
//   d=ndim f=input_file e=nev p=npr q=npr1 t=elemTypeStr w=targ_en
//   x=nx y=ny z=nz 0=x0 1=x1 2=y0 3=y1 4=z0 5=z1 6=pbc_x 7=pbc_y 8=pbc_z
  while ((opt = getopt(argc, argv, "f:d:e:p:q:t:x:y:z:0:1:2:3:4:5:6:7:8:")) != -1) {
    switch (opt) {
      case 'f':
        inputFile = optarg;
        break;
      case 'd':
        variables["ndim"] = optarg;
        break;
      case 'e':
        variables["nev"] = optarg;
        break;
      case 'p':
        variables["npr"] = optarg;
        break;
      case 'q':
        variables["npr1"] = optarg;
        break;
      case 't':
        variables["elemTypeStr"] = optarg;
        break;
      case 'w':
        variables["targ_en"] = optarg;
        break;
      case 'x':
        variables["nx"] = optarg;
        break;
      case 'y':
        variables["ny"] = optarg;
        break;
      case 'z':
        variables["nz"] = optarg;
        break;
      case '0':
        variables["x0"] = optarg;
        break;
      case '1':
        variables["x1"] = optarg;
        break;
      case '2':
        variables["y0"] = optarg;
        break;
      case '3':
        variables["y1"] = optarg;
        break;
      case '4':
        variables["z0"] = optarg;
        break;
      case '5':
        variables["z1"] = optarg;
        break;
      case '6':
        variables["pbc_x"] = optarg;
        break;
      case '7':
        variables["pbc_y"] = optarg;
        break;
      case '8':
        variables["pbc_z"] = optarg;
        break;
      default:
        std::cerr << "Usage: " << argv[0] << " [-f input_file] [-d ndim] [-e nev] [-p npr] [-q npr1] [t elemTypeStr] [-w targ_en] [-x nx] [-y ny] [-z nz] [-0 x0] [-1 x1] [-2 y0] [-3 y1] [-4 z0] [-5 z1] [-6 pbc_x] [-7 pbc_y] [-8 pbc_z]" << std::endl;
        return 1;
    }
  }

  if (!inputFile.empty()) {
    parseFile(inputFile);
  }

// wasteful, but assign values of the variables array to the proper variable names
ndim=stoi(variables["ndim"]);
nev=stoi(variables["nev"]);
npr=stoi(variables["npr"]);
nx=stoi(variables["nx"]);
targ_en=stod(variables["targ_en"]);
x0=stod(variables["x0"]);
x1=stod(variables["x1"]);
pbc_x=stoi(variables["pbc_x"]);
ny=stoi(variables["ny"]);
y0=stod(variables["y0"]);
y1=stod(variables["y1"]);
pbc_y=stoi(variables["pbc_y"]);
nz=stoi(variables["nz"]);
z0=stod(variables["z0"]);
z1=stod(variables["z1"]);
elemTypeStr=variables["elemTypeStr"];
npr1=stoi(variables["npr1"]);
pbc_z=stoi(variables["pbc_z"]);

  // Because libMesh::ElemType is an enumeration type, I can't directly assign a string value to it.
  // Instead we use this mapping to convert an input strong elemTypeStr to the enum libMesh::ElemType
  std::unordered_map<std::string, libMesh::ElemType> elemTypeMap = {
    {"EDGE2", libMesh::EDGE2},  {"EDGE3", libMesh::EDGE3},  {"EDGE4", libMesh::EDGE4},
    {"TRI3", libMesh::TRI3},  {"TRI6", libMesh::TRI6},
    {"QUAD4", libMesh::QUAD4},  {"QUAD8", libMesh::QUAD8},  {"QUAD9", libMesh::QUAD9},
    {"TET4", libMesh::TET4},{"TET10", libMesh::TET10},
    {"HEX8", libMesh::HEX8}, {"HEX20", libMesh::HEX20}, {"HEX27", libMesh::HEX27},
    {"PRISM6", libMesh::PRISM6},  {"PRISM15", libMesh::PRISM15},  {"PRISM18", libMesh::PRISM18},
    {"PYRAMID5", libMesh::PYRAMID5},  {"PYRAMID13", libMesh::PYRAMID13},  {"PYRAMID14", libMesh::PYRAMID14}
  };
  // Convert the string to libMesh::ElemType
    auto it = elemTypeMap.find(elemTypeStr);
    if (it == elemTypeMap.end())
    {
        std::cerr << "Invalid element type: " << elemTypeStr << std::endl;
        return 1;
    }
    libMesh::ElemType elemType = it->second;

  if (ndim == 1)
    logfile << "elemTypeStr, ndim, nev, nx = " << elemTypeStr << " " << ndim << " " << nev << " " << nx << "\n";
    std::cout << "elemTypeStr, ndim, nev, nx = " << elemTypeStr << " " << ndim << " " << nev << " " << nx << "\n";
  if (ndim == 2)
    logfile << "elemTypeStr, ndim, nev, nx, ny = " << elemTypeStr << " " << ndim << " " << nev << " " << nx << " " << ny << "\n";
    std::cout << "elemTypeStr, ndim, nev, nx, ny = " << elemTypeStr << " " << ndim << " " << nev << " " << nx << " " << ny << "\n";
  if (ndim == 3)
    logfile << "elemTypeStr, ndim, nev, nx, ny, nz = " << elemTypeStr << " " << ndim << " " << nev << " " << nx << " " << ny << " " << nz << "\n";
    std::cout << "elemTypeStr, ndim, nev, nx, ny, nz = " << elemTypeStr << " " << ndim << " " << nev << " " << nx << " " << ny << " " << nz << "\n";

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

//  libmesh_example_requires(false, "--enable-periodic");

  // Create a mesh, with dimension to be overridden later, on the
  // default MPI communicator.
  Mesh mesh(init.comm());

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);
	MeshRefinement mesh_refinement (mesh);

  // Create a EigenSystem named "Eigensystem" and (for convenience)
  // use a reference to the system we create.
  EigenSystem & eigen_system =
    equation_systems.add_system<EigenSystem> ("Eigensystem");

        logfile << "periodic boundary conditions flags:" << pbc_x << " " 
	          << pbc_y << " " << pbc_z << std::endl;
        std::cout << "periodic boundary conditions flags:" << pbc_x << " " 
	          << pbc_y << " " << pbc_z << std::endl;
  // PERIODIC BOUNDARY CONDITIONS 
    DofMap & dof_map = eigen_system.get_dof_map();
    if ((ndim == 1) && (pbc_x == 1)){
        logfile << "periodic boundary conditions for x" << std::endl;
        std::cout << "periodic boundary conditions for x" << std::endl;
	PeriodicBoundary xbdry(RealVectorValue(x1-x0));
        xbdry.myboundary = 0;
        xbdry.pairedboundary = 1;
        eigen_system.get_dof_map().add_periodic_boundary(xbdry);
    }

    if ((ndim == 2) && ((pbc_x == 1) || (pbc_y == 1))){
        logfile << "periodic boundary conditions for " ;
        std::cout << "periodic boundary conditions for " ;

	if ((pbc_x == 1) && (pbc_y == 1)){
          logfile << "x and y" << std::endl;
          std::cout << "x and y" << std::endl;

          PeriodicBoundary xbdry(RealVectorValue(x1-x0, 0.0)); 
          xbdry.myboundary = 3;
          xbdry.pairedboundary = 1;
          eigen_system.get_dof_map().add_periodic_boundary(xbdry); 

          PeriodicBoundary ybdry(RealVectorValue(0.0, y1-y0));
	  ybdry.myboundary = 0;
          ybdry.pairedboundary = 2;
          eigen_system.get_dof_map().add_periodic_boundary(ybdry); 
	}
	if ((pbc_x == 0) && (pbc_y == 1)){
          logfile << "y" << std::endl;
          std::cout << "y" << std::endl;

          PeriodicBoundary xbdry(RealVectorValue(x1-x0, 0.0)); 
          xbdry.myboundary = 3;
          xbdry.pairedboundary = 1;
          eigen_system.get_dof_map().add_periodic_boundary(xbdry); 
	}
	if ((pbc_x == 1) && (pbc_y == 0)){
          logfile << "x" << std::endl;
          std::cout << "x" << std::endl;

          PeriodicBoundary ybdry(RealVectorValue(0.0, y1-y0));
	  ybdry.myboundary = 0;
          ybdry.pairedboundary = 2;
          eigen_system.get_dof_map().add_periodic_boundary(ybdry); 
	}
    }

    if ((ndim == 3) && ((pbc_x == 1) || (pbc_y == 1) || (pbc_z == 1))){
        logfile << "periodic boundary conditions for " ;
        std::cout << "periodic boundary conditions for " ;

	if ((pbc_x == 1) && (pbc_y == 1) && (pbc_z == 1)){
          logfile << "x, y, and z " << std::endl;
          std::cout << "x, y, and z " << std::endl;

          PeriodicBoundary xbdry(RealVectorValue(x1-x0, 0.0, 0.0)); 
          xbdry.myboundary = 4;
          xbdry.pairedboundary = 2;
          eigen_system.get_dof_map().add_periodic_boundary(xbdry); 

          PeriodicBoundary ybdry(RealVectorValue(0.0, y1-y0, 0.0));
	  ybdry.myboundary = 1;
          ybdry.pairedboundary = 3;
          eigen_system.get_dof_map().add_periodic_boundary(ybdry); 

          PeriodicBoundary zbdry(RealVectorValue(0.0, 0.0, z1-z0));
	  zbdry.myboundary = 0;
          zbdry.pairedboundary = 5;
          eigen_system.get_dof_map().add_periodic_boundary(zbdry); 
	}
	if ((pbc_x == 1) && (pbc_y == 1) && (pbc_z == 0)){
          logfile << "x and y " << std::endl;
          std::cout << "x and y " << std::endl;

          PeriodicBoundary xbdry(RealVectorValue(x1-x0, 0.0, 0.0)); 
          xbdry.myboundary = 4;
          xbdry.pairedboundary = 2;
          eigen_system.get_dof_map().add_periodic_boundary(xbdry); 

          PeriodicBoundary ybdry(RealVectorValue(0.0, y1-y0, 0.0));
	  ybdry.myboundary = 1;
          ybdry.pairedboundary = 3;
          eigen_system.get_dof_map().add_periodic_boundary(ybdry); 
	}
	if ((pbc_x == 1) && (pbc_y == 0) && (pbc_z == 1)){
          logfile << "x and z " << std::endl;
          std::cout << "x and z " << std::endl;

          PeriodicBoundary xbdry(RealVectorValue(x1-x0, 0.0, 0.0)); 
          xbdry.myboundary = 4;
          xbdry.pairedboundary = 2;
          eigen_system.get_dof_map().add_periodic_boundary(xbdry); 

          PeriodicBoundary zbdry(RealVectorValue(0.0, 0.0, z1-z0));
	  zbdry.myboundary = 0;
          zbdry.pairedboundary = 5;
          eigen_system.get_dof_map().add_periodic_boundary(zbdry); 
	}
	if ((pbc_x == 0) && (pbc_y == 1) && (pbc_z == 1)){
          logfile << "y and z " << std::endl;
          std::cout << "y and z " << std::endl;

          PeriodicBoundary ybdry(RealVectorValue(0.0, y1-y0, 0.0));
	  ybdry.myboundary = 1;
          ybdry.pairedboundary = 3;
          eigen_system.get_dof_map().add_periodic_boundary(ybdry); 

          PeriodicBoundary zbdry(RealVectorValue(0.0, 0.0, z1-z0));
	  zbdry.myboundary = 0;
          zbdry.pairedboundary = 5;
          eigen_system.get_dof_map().add_periodic_boundary(zbdry); 
	}
	if ((pbc_x == 1) && (pbc_y == 0) && (pbc_z == 0)){
          logfile << "x" << std::endl;
          std::cout << "x" << std::endl;

          PeriodicBoundary xbdry(RealVectorValue(x1-x0, 0.0, 0.0)); 
          xbdry.myboundary = 4;
          xbdry.pairedboundary = 2;
          eigen_system.get_dof_map().add_periodic_boundary(xbdry); 
	}
	if ((pbc_x == 0) && (pbc_y == 1) && (pbc_z == 0)){
          logfile << "y" << std::endl;
          std::cout << "y" << std::endl;

          PeriodicBoundary ybdry(RealVectorValue(0.0, y1-y0, 0.0));
	  ybdry.myboundary = 1;
          ybdry.pairedboundary = 3;
          eigen_system.get_dof_map().add_periodic_boundary(ybdry); 
	}
	if ((pbc_x == 0) && (pbc_y == 0) && (pbc_z == 1)){
          logfile << "z" << std::endl;
          std::cout << "z" << std::endl;

          PeriodicBoundary zbdry(RealVectorValue(0.0, 0.0, z1-z0));
	  zbdry.myboundary = 0;
          zbdry.pairedboundary = 5;
          eigen_system.get_dof_map().add_periodic_boundary(zbdry); 
	}
    }

  // Use the internal mesh generator to create a uniform
  // 1D grid on a line.
  if (ndim == 1)
  MeshTools::Generation::build_line (mesh, nx, x0, x1,
                                       elemType);       // EDGE2 EDGE3 EDGE4
  // 2D grid on a square.
  else if (ndim == 2)
  MeshTools::Generation::build_square (mesh, nx, ny, x0, x1, y0, y1,
                                       elemType);       // TRI3 TRI6 QUAD4 QUAD8 QUAD9
  // 3D grid on a cube.
  else if (ndim == 3)
  MeshTools::Generation::build_cube (mesh, nx, ny, nz, x0, x1, y0, y1, z0, z1,
                                       elemType);    // TET4 TET10 HEX8 HEX20 HEX27 PRISM6 PRISM15 PRISM18 PYRAMID5 PYRAMID13 PYRAMID14
  // ADD TRAP IF ndim != 1,2,3

  // Print information about the mesh to the screen.
  mesh.print_info(logfile);
  mesh.print_info(std::cout);

        mesh_refinement.set_periodic_boundaries_ptr(dof_map.get_periodic_boundaries());


int n_nodes = mesh.n_nodes();
int n_elem = mesh.n_elem();
logfile << "# nodes = " << n_nodes << "\n";
logfile << "# elements = " << n_elem << "\n";
logfile << "Writing nodal xyz coordinates to file: results/femvib.xyz" << std::endl;
std::cout << "# nodes = " << n_nodes << "\n";
std::cout << "# elements = " << n_elem << "\n";
std::cout << "Writing nodal xyz coordinates to file: results/femvib.xyz" << std::endl;
std::ofstream xyzfile;
xyzfile.open("results/femvib.xyz");

MeshBase::const_node_iterator node_it = mesh.nodes_begin();
const MeshBase::const_node_iterator node_end = mesh.nodes_end();

xyzfile << "Number of nodes: " << mesh.n_nodes() << "\n";

int i_node {0};
for (; node_it != node_end; ++node_it)
{
    const libMesh::Point& p = **node_it;
    const libMesh::Node* node = *node_it;
//    const libMesh::Point& p = node->point();
//
    const Real px = p(0);
    if (ndim == 1){ 
      const Real pw = p(1);
      xyzfile << px << "\n";
    }
    else if (ndim == 2){
      const Real py = p(1);
      const Real pw = p(2);
      xyzfile << px << " " << py << "\n";
    }
    else if (ndim == 3){
      const Real py = p(1);
      const Real pz = p(2);
      const Real pw = p(3);
      xyzfile << px << " " << py << " " << pz << "\n";
    }
	++i_node;
//	std::cout << i_node << std::endl;
}

    xyzfile.close();
////    std::cout << "Finished writing nodal xyz coordinates to file: " << "%22%3E%3Cimg%20src.xyz" << std::endl;
//    std::cout << "Finished writing nodal xyz coordinates to file: " << "results/femvib.xyz" << std::endl;
//

logfile << "Writing number of dimensions to file: results/ndim" << std::endl;
std::cout << "Writing number of dimensions to file: results/ndim" << std::endl;
  std::ofstream ndimfile;
  ndimfile.open("results/ndim");
  ndimfile << ndim ;
  ndimfile.close();

  // Declare the system variables.
  // Adds the variable "p" to "Eigensystem".   "p"
  // will be approximated using second-order approximation.
  eigen_system.add_variable("p", FIRST);
  eigen_system.attach_assemble_function (assemble_mass);
  equation_systems.parameters.set<unsigned int>("eigenpairs")    = nev;
  equation_systems.parameters.set<unsigned int>("basis vectors") = nev*3;
      eigen_system.eigen_solver->set_eigensolver_type(KRYLOVSCHUR);
      eigen_system.eigen_solver->set_position_of_spectrum(SMALLEST_REAL);
  equation_systems.parameters.set<Real>("linear solver tolerance") = pow(TOLERANCE, 5./3.);
  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 1000;
  eigen_system.set_eigenproblem_type(GHEP);
  equation_systems.init();
  eigen_system.assemble();

  eigen_system.get_matrix_A().close();  // Finalize stiffness matrix
  eigen_system.get_matrix_B().close();  // Finalize mass matrix

// Downcast matrix wrappers to access PETSc raw matrices
auto& A_mat_wrapper = static_cast<libMesh::PetscMatrix<Number>&>(eigen_system.get_matrix_A());
auto& B_mat_wrapper = static_cast<libMesh::PetscMatrix<Number>&>(eigen_system.get_matrix_B());

// Extract raw PETSc matrix handles
Mat A = A_mat_wrapper.mat();
Mat B = B_mat_wrapper.mat();

  equation_systems.print_info(std::cout);
  equation_systems.print_info(logfile);

  EPS eps;
EPSCreate(PETSC_COMM_WORLD, &eps);

EPSSetOperators(eps, A, B);
EPSSetProblemType(eps, EPS_GHEP);  // Generalized Hermitian

EPSSetType(eps, EPSKRYLOVSCHUR);  // Or any other appropriate solver
EPSSetWhichEigenpairs(eps, EPS_TARGET_REAL);
EPSSetTarget(eps, targ_en);  // Your desired energy target

// Enable shift-and-invert
ST st;
EPSGetST(eps, &st);
STSetType(st, STSINVERT);

// Set dimensions (nev = how many eigenpairs to find)
EPSSetDimensions(eps, nev, PETSC_DEFAULT, PETSC_DEFAULT);

// Finalize and solve
EPSSetFromOptions(eps);
EPSSolve(eps);

// Optional: inspect spectrum
PetscInt nconv;
EPSGetConverged(eps, &nconv);

// Problem: When using target_energy, SLEPc eigenpairs are unsorted. 
// We will sort on the real part of the eigenvalue and asume the imaginary part = 0.
// Extract all eigenvalues and their original indices
std::vector<std::pair<double, PetscInt>> eigenpairs;

for (PetscInt i = 0; i < nconv; ++i)
{
    PetscScalar kr, ki;
    EPSGetEigenvalue(eps, i, &kr, &ki);

    if (std::abs(PetscImaginaryPart(ki)) > 1e-10) continue;  // skip complex eigenvalues

    eigenpairs.emplace_back(PetscRealPart(kr), i);
}

// Sort by increasing eigenvalue
std::sort(eigenpairs.begin(), eigenpairs.end(),
          [](const auto& a, const auto& b) {
              return a.first < b.first;
          });

// Loop over sorted eigenvalues 
std::ofstream eigenvalue_file("results/eigenvalues.txt");
eigenvalue_file << "\n EIGENVALUES\ni   Re(evalue) \n";
eigenvalue_file << std::fixed << std::setprecision(5);
logfile << "\n EIGENVALUES\ni    Re(evalue)    Eigenstate file\n";
logfile << std::fixed << std::setprecision(5);
for (std::size_t i = 0; i < eigenpairs.size(); ++i)
{
    double eval = eigenpairs[i].first;
    PetscInt orig_index = eigenpairs[i].second;

    // print eigenstate if npr1 <= i < npr1+npr
    std::ostringstream eigenstate_filename;
    if (npr1 <= i && i < npr1 + npr) {
      eigenstate_filename << "results/eigenstate_" << i << ".txt";
      write_petsc_eigenvector_ascii(eps, orig_index, equation_systems,
                                   "Eigensystem", "p", eigenstate_filename.str());
    }
    // print eigenvalue (real part only)
    eigenvalue_file << i << "  " << eval << std::endl;
    logfile << i << "  " << eval << "  " << eigenstate_filename.str() << std::endl;
}
  std::cout << "Wrote " << npr << " eigenstates to results/\n" << std::endl;

  // All done.
  auto end = std::chrono::system_clock::now();  // end timekeeping

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);

  logfile << "\n finished computation at " << std::ctime(&end_time)
              << std::fixed << std::setprecision(3)
              << "Elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
  std::cout << "finished computation at " << std::ctime(&end_time)
              << std::fixed << std::setprecision(3)
              << "Elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;

return 0;
}

void assemble_mass(EquationSystems & es,
                   const std::string & libmesh_dbg_var(system_name))
{
 using namespace std;
 int i;
// necessary to initialize nptspot below
 int nptspot {0}, nptsx, nptsxx, nptsxy, nptsyy, nptszx, nptszy, nptszz;

 double avpot;
 double sdpot;
 double adev, varia, skew, kurt;
 int subtracpot;
 int sonpot;
 int ndim;
 double cut=0.0, cutxx, cutxy, cutyy, cutzx, cutzy, cutzz, rval;

// std::ifstream inpot;
 std::ifstream inpot, inxx, inxy, inyy, inzx, inzy, inzz;
 std::ofstream outpot;
 std::ofstream logfile("results/femvib.log", std::ios::app);

logfile << "Reading Kriging parameters from files results/optr results/opp" << std::endl;
std::cout << "Reading Kriging parameters from files results/optr results/opp" << std::endl;
  std::ifstream oppmfile;
  oppmfile.open("results/opp");
  oppmfile >> nptsx;
  oppmfile.close();
  std::ifstream optrfile;
  optrfile.open("results/optr");
  optrfile >> rval;
  optrfile.close();
 logfile << "rval = " << rval << "\n" ;
 std::cout << "rval = " << rval << "\n" ;
 logfile << "nptsx = " << nptsx << "\n" ;
 std::cout << "nptsx = " << nptsx << "\n" ;

logfile << "Reading number of dimensions from file results/ndim" << std::endl;
std::cout << "Reading number of dimensions from file results/ndim" << std::endl;
  std::ifstream ndimfile;
  ndimfile.open("results/ndim");
  ndimfile >> ndim ;
  ndimfile.close();

// DEBUG 
// fstream outAm;
// outAm.open("results/A_m.out", ios::out);
// outAm << fixed;	// fixes the precision of the output including trailing zeros
// outAm.precision(3);	// fixes the precision of the output to 3 decimal places
// DEBUG 

 inpot.open("results/energyP");
 outpot.open("results/pot.out");
 outpot << fixed;	// fixes the precision of the output including trailing zeros
 outpot.precision(3);	// fixes the precision of the output to 3 decimal places
 // The precision can be a conditional in fem.sh determined by the input data
 if (!inpot)
 {
	std::cerr << "Error: file energyP could not be opened\n";
	exit(1);
 }
 std::string number;
 std::string dataline;
 inpot >> nptspot;  // 
 logfile << "Number of energies in PES grid = " << nptspot << "\n";
 std::cout << "Number of energies in PES grid = " << nptspot << "\n";

 Eigen::MatrixXd m(nptspot,ndim);
 Eigen::VectorXd v(nptspot);

 MatDoub xpot(nptspot,ndim);
 VecDoub Energy(nptspot);
 
 for (i=0;i<nptspot;i++)
 {
	cout.precision(10);
	cout.setf(ios::fixed,ios::floatfield);
if (ndim == 1)	inpot >> xpot[i][0] >> Energy[i];
if (ndim == 2)	inpot >> xpot[i][0] >> xpot[i][1] >> Energy[i];
if (ndim == 3)	inpot >> xpot[i][0] >> xpot[i][1] >> xpot[i][2] >> Energy[i];
 }
 logfile << "\n" ;
 std::cout << "\n" ;

 outpot.close();
 inpot.close();

 std::unique_ptr<Shep_interp> krigxx_ptr = nullptr;
 std::unique_ptr<Shep_interp> krigxy_ptr = nullptr;
 std::unique_ptr<Shep_interp> krigyy_ptr = nullptr;
 std::unique_ptr<Shep_interp> krigzx_ptr = nullptr;
 std::unique_ptr<Shep_interp> krigzy_ptr = nullptr;
 std::unique_ptr<Shep_interp> krigzz_ptr = nullptr;
 std::unique_ptr<MatDoub> xxx_ptr = nullptr;
 std::unique_ptr<MatDoub> xxy_ptr = nullptr;
 std::unique_ptr<MatDoub> xyy_ptr = nullptr;
 std::unique_ptr<MatDoub> xzx_ptr = nullptr;
 std::unique_ptr<MatDoub> xzy_ptr = nullptr;
 std::unique_ptr<MatDoub> xzz_ptr = nullptr;
 std::unique_ptr<VecDoub> xxvals_ptr = nullptr;
 std::unique_ptr<VecDoub> xyvals_ptr = nullptr;
 std::unique_ptr<VecDoub> yyvals_ptr = nullptr;
 std::unique_ptr<VecDoub> zxvals_ptr = nullptr;
 std::unique_ptr<VecDoub> zyvals_ptr = nullptr;
 std::unique_ptr<VecDoub> zzvals_ptr = nullptr;
 
// read G_xx matrix elements
 inxx.open("results/matrixg.xx");
                 if (!inxx)
                        {
                                std::cerr << "Error: file xx could not be opened\n";
                                exit(1);
                        }
 inxx >> nptsxx;
// MatDoub xxx(nptsxx,ndim);
// VecDoub xxvals(nptsxx); 
 xxx_ptr = std::make_unique<MatDoub>(nptsxx,ndim);
 xxvals_ptr = std::make_unique<VecDoub>(nptsxx); 

 for (i=0;i<nptsxx;i++)
 {
	cout.precision(10);
	cout.setf(ios::fixed,ios::floatfield);
	if (ndim == 1)	inxx >> (*xxx_ptr)[i][0] >> (*xxvals_ptr)[i];
	if (ndim == 2)	inxx >> (*xxx_ptr)[i][0] >> (*xxx_ptr)[i][1] >> (*xxvals_ptr)[i];
	if (ndim == 3)	inxx >> (*xxx_ptr)[i][0] >> (*xxx_ptr)[i][1] >> (*xxx_ptr)[i][2] >> (*xxvals_ptr)[i];
 }
 inxx.close();

 if (ndim > 1){
// read G_xy matrix elements
 inxy.open("results/matrixg.xy");
                 if (!inxy)
                        {
                                std::cerr << "Error: file xy could not be opened\n";
                                exit(1);
                        }
 inxy >> nptsxy;
// MatDoub xxy(nptsxy,ndim);
// VecDoub xyvals(nptsxy);
 xxy_ptr = std::make_unique<MatDoub>(nptsxy,ndim);
 xyvals_ptr = std::make_unique<VecDoub>(nptsxy); 
 for (i=0;i<nptsxy;i++)
 {
	cout.precision(10);
	cout.setf(ios::fixed,ios::floatfield);
	if (ndim == 2)	inxy >> (*xxy_ptr)[i][0] >> (*xxy_ptr)[i][1] >> (*xyvals_ptr)[i];
	if (ndim == 3)	inxy >> (*xxy_ptr)[i][0] >> (*xxy_ptr)[i][1] >> (*xxy_ptr)[i][2] >> (*xyvals_ptr)[i];
 }
 inxy.close();

// read G_yy matrix elements
 inyy.open("results/matrixg.yy");
                 if (!inyy)
                        {
                                std::cerr << "Error: file yy could not be opened\n";
                                exit(1);
                        }
 inyy >> nptsyy;
// MatDoub xyy(nptsyy,ndim);
// VecDoub yyvals(nptsyy);
 xyy_ptr = std::make_unique<MatDoub>(nptsyy,ndim);
 yyvals_ptr = std::make_unique<VecDoub>(nptsyy); 
 for (i=0;i<nptsyy;i++)
 {
	cout.precision(10);
	cout.setf(ios::fixed,ios::floatfield);
	if (ndim == 2)	inyy >> (*xyy_ptr)[i][0] >> (*xyy_ptr)[i][1] >> (*yyvals_ptr)[i];
	if (ndim == 3)	inyy >> (*xyy_ptr)[i][0] >> (*xyy_ptr)[i][1] >> (*xyy_ptr)[i][2] >> (*yyvals_ptr)[i];
 }
 inyy.close();
 }

 if (ndim > 2){
// read G_zx matrix elements
 inzx.open("results/matrixg.zx");
                 if (!inzx)
                        {
                                std::cerr << "Error: file zx could not be opened\n";
                                exit(1);
                        }
 inzx >> nptszx;
// MatDoub xzx(nptszx,ndim);
// VecDoub zxvals(nptszx);
 xzx_ptr = std::make_unique<MatDoub>(nptszx,ndim);
 zxvals_ptr = std::make_unique<VecDoub>(nptszx); 
 for (i=0;i<nptszx;i++)
 {
	cout.precision(10);
	cout.setf(ios::fixed,ios::floatfield);
	inzx >> (*xzx_ptr)[i][0] >> (*xzx_ptr)[i][1] >> (*xzx_ptr)[i][2] >> (*zxvals_ptr)[i];
 }
 inzx.close();

// read G_zy matrix elements
 inzy.open("results/matrixg.zy");
                 if (!inzy)
                        {
                                std::cerr << "Error: file zy could not be opened\n";
                                exit(1);
                        }
 inzy >> nptszy;
// MatDoub xzy(nptszy,ndim);
// VecDoub zyvals(nptszy);
 xzy_ptr = std::make_unique<MatDoub>(nptszy,ndim);
 zyvals_ptr = std::make_unique<VecDoub>(nptszy); 
 for (i=0;i<nptszy;i++)
 {
	cout.precision(10);
	cout.setf(ios::fixed,ios::floatfield);
	inzy >> (*xzy_ptr)[i][0] >> (*xzy_ptr)[i][1] >> (*xzy_ptr)[i][2] >> (*zyvals_ptr)[i];
 }
 inzy.close();

// read G_zz matrix elements
 inzz.open("results/matrixg.zz");
                 if (!inzz)
                        {
                                std::cerr << "Error: file zz could not be opened\n";
                                exit(1);
                        }
 inzz >> nptszz;
// MatDoub xzz(nptszz,ndim);
// VecDoub zzvals(nptszz);
 xzz_ptr = std::make_unique<MatDoub>(nptszz,ndim);
 zzvals_ptr = std::make_unique<VecDoub>(nptszz); 
 for (i=0;i<nptszz;i++)
 {
	cout.precision(10);
	cout.setf(ios::fixed,ios::floatfield);
	inzz >> (*xzz_ptr)[i][0] >> (*xzz_ptr)[i][1] >> (*xzz_ptr)[i][2] >> (*zzvals_ptr)[i];
 }
 inzz.close();
 }

 // Create the Shep_interp objects; "krig??" are being created here as instances of the Shep_interp struct.
//  Shep_interp krigxx(xxx,xxvals);
//  Shep_interp krigxy(xxy,xyvals);
//  Shep_interp krigyy(xyy,yyvals);
//  Shep_interp krigzx(xzx,zxvals);
//  Shep_interp krigzy(xzy,zyvals);
//  Shep_interp krigzz(xzz,zzvals);
   krigxx_ptr = std::make_unique<Shep_interp>((*xxx_ptr), (*xxvals_ptr));
   if (ndim > 1){
     krigyy_ptr = std::make_unique<Shep_interp>((*xyy_ptr), (*yyvals_ptr));
     krigxy_ptr = std::make_unique<Shep_interp>((*xxy_ptr), (*xyvals_ptr));
   }
   if (ndim > 2){
     krigzx_ptr = std::make_unique<Shep_interp>((*xzx_ptr), (*zxvals_ptr));
     krigzy_ptr = std::make_unique<Shep_interp>((*xzy_ptr), (*zyvals_ptr));
     krigzz_ptr = std::make_unique<Shep_interp>((*xzz_ptr), (*zzvals_ptr));
   }

  fstream potfile;
  potfile.open("results/pot.out", ios::out);
  if(!potfile) {cerr << "Error opening PES output file."; return;}

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Eigensystem");

#ifdef LIBMESH_HAVE_SLEPC

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running.
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to our system.
  EigenSystem & eigen_system = es.get_system<EigenSystem> ("Eigensystem");

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = eigen_system.get_dof_map().variable_type(0);

  SparseMatrix<Number> & matrix_A = eigen_system.get_matrix_A();
  SparseMatrix<Number> & matrix_B = eigen_system.get_matrix_B();


  // A reference to the two system matrices
  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a std::unique_ptr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));

  // A  Gauss quadrature rule for numerical integration.
  // Use the default quadrature order.
  QGauss qrule (dim, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // from FEMvib:
  std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));
  QGauss qface (dim-1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule (&qface);

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // ALC -- added from 140522_19_07_51.C
  const std::vector<Point>& q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real>> & phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const DofMap & dof_map = eigen_system.get_dof_map();

  // The element mass and stiffness matrices.
  DenseMatrix<Number> Me;
  DenseMatrix<Number> Ke;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  MeshBase::const_element_iterator       el     = mesh.elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.elements_end();
  logfile << "Writing PES data to file: results/pot.out " << std::endl;
  std::cout << "Writing PES data to file: results/pot.out " << std::endl;

// Check for duplicate points in xpot -- these will wreck the numerical integration
  dupl(xpot,nptspot,ndim);

  double V, gxx, gxy, gyy, gzx, gzy, gzz;
  MatDoub pts(nptsx,ndim);
  VecDoub yvecnew(nptsx);


  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  In case users
  // later modify this program to include refinement, we will
  // be safe and will only consider the active elements;
  // hence we use a variant of the active_elem_iterator.

//  for (const auto & elem : mesh.active_local_element_ptr_range())
  for ( ; el != end_el; ++el)   // is this the same as the above line?  Above is from eigenproblem example.
    {
      const Elem* elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);

      // Zero the element matrices before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      const unsigned int n_dofs =
        cast_int<unsigned int>(dof_indices.size());
      Ke.resize (n_dofs, n_dofs);
      Me.resize (n_dofs, n_dofs);

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      //
      // We will build the element matrix.  This involves
      // a double loop to integrate the test functions (i) against
      // the trial functions (j).
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	{
	const double x = q_point[qp](0);
        const double y = q_point[qp](1);
        const double z = q_point[qp](2);
        VecDoub potpoint(ndim);
        potpoint[0] =  x;
        potpoint[1] =  y;
        potpoint[2] =  z;
	// there's a function "howmany" defined elsewhere that takes 2 args
        howmany8(xpot, Energy, potpoint, ndim, nptspot, nptsx, pts, yvecnew);

        Powvargram vgram(pts,yvecnew,rval);
        Krig<Powvargram> krigp(pts,yvecnew,vgram);

        V = krigp.interp(potpoint);

//        gxx = krigxx.interp(potpoint);
//        gxy = krigxy.interp(potpoint);
//        gyy = krigyy.interp(potpoint);
//        gzx = krigzx.interp(potpoint);
//        gzy = krigzy.interp(potpoint);
//        gzz = krigzz.interp(potpoint);
           gxx = krigxx_ptr->interp(potpoint);
	if (ndim > 1){
           gxy = krigxy_ptr->interp(potpoint);
           gyy = krigyy_ptr->interp(potpoint);
        }
	if (ndim > 2){
           gzx = krigzx_ptr->interp(potpoint);
           gzy = krigzy_ptr->interp(potpoint);
           gzz = krigzz_ptr->interp(potpoint);
	}

                        cout.precision(10);
                        cout.setf(ios::fixed,ios::floatfield);
        if (ndim == 1) potfile << x << '\t' << V << '\t' <<endl;
        if (ndim == 2) potfile << x << '\t' << y << '\t' << V << '\t' <<endl;
        if (ndim == 3) potfile << x << '\t' << y << '\t' << z <<  '\t' << V <<'\t' <<endl; 
			//
        for (unsigned int i=0; i<n_dofs; i++)
          for (unsigned int j=0; j<n_dofs; j++)
            {
              Me(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
//              Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);  // this line from eigenproblems example
              Ke(i,j) += JxW[qp]*(
                              0.5*(   gxx*dphi[i][qp](0)*dphi[j][qp](0)
                                  )*33.71526 + V*phi[i][qp]*phi[j][qp]);
	      if (ndim > 1){
                Ke(i,j) += JxW[qp]*(
                              0.5*(     gxy*dphi[i][qp](0)*dphi[j][qp](1)
                                      + gxy*dphi[i][qp](1)*dphi[j][qp](0)
                                      + gyy*dphi[i][qp](1)*dphi[j][qp](1)
                                  )*33.71526);
	      }
	      if (ndim > 2){
                Ke(i,j) += JxW[qp]*(
                              0.5*(     gzz*dphi[i][qp](2)*dphi[j][qp](2)
                                      + gzx*dphi[i][qp](0)*dphi[j][qp](2)
                                      + gzx*dphi[i][qp](2)*dphi[j][qp](0)
                                      + gzy*dphi[i][qp](1)*dphi[j][qp](2)
                                      + gzy*dphi[i][qp](2)*dphi[j][qp](1)
                                  )*33.71526);
	      }
                                   //   + dgxxdx*dphi[i][qp](0)
                                   //   + dgyydy*dphi[i][qp](1)
                                   //   + dgzzdz*dphi[i][qp](2)
                                   //   + dgxydx*dphi[i][qp](1)
                                   //   + dgxydy*dphi[i][qp](0)
                                   //   + dgxzdx*dphi[i][qp](2)
                                   //   + dgzxdz*dphi[i][qp](0)
                                   //   + dgyzdy*dphi[i][qp](2)
                                  //    + dgzydz*dphi[i][qp](1)
//                                  )*33.71526 + V*phi[i][qp]*phi[j][qp]);
//  DEBUG	
//	      outAm << i << '\t' << j << '\t' << Ke(i,j) << '\t' << phi[i][qp] << '\t' 
//	            << phi[j][qp] << '\t' << V << '\t' << gxx << '\t' << gxy << '\t' 
//		    << gyy << endl;
//  DEBUG	
            }
        }

      // On an unrefined mesh, constrain_element_matrix does
      // nothing.  If this assembly function is ever repurposed to
      // run on a refined mesh, getting the hanging node constraints
      // right will be important.  Note that, even with
      // asymmetric_constraint_rows = false, the constrained dof
      // diagonals still exist in the matrix, with diagonal entries
      // that are there to ensure non-singular matrices for linear
      // solves but which would generate positive non-physical
      // eigenvalues for eigensolves.
      // THIS STEP IS NECESSARY FOR PERIODIC BOUNDARY CONDITIONS
       dof_map.constrain_element_matrix(Ke, dof_indices, false);
       dof_map.constrain_element_matrix(Me, dof_indices, false);

      // Finally, simply add the element contribution to the
      // overall matrices A and B.
      matrix_A.add_matrix (Ke, dof_indices);
      matrix_B.add_matrix (Me, dof_indices);
    } // end of element loop

      potfile.close();
// DEBUG      outAm.close();
      logfile << "Finished writing PES data to file: " << "results/pot.out" << std::endl;
      std::cout << "Finished writing PES data to file: " << "results/pot.out" << std::endl;
      logfile << "Matrix A and B assembly completed. " << "\n" << std::endl;
      std::cout << "Matrix A and B assembly completed. " << "\n" << std::endl;

      //      THE FOLLOWING COMMANDS GENERATE A DEPRECATED CODE WARNING
      //      but still (as of petsc-3.23.4) geneate A.m and B.m files if needed for debugging
      //      These are big files so better left out if not needed
//      matrix_A.print_matlab("results/A.m");
//      matrix_B.print_matlab("results/B.m");
//      std::cout << "Matrix A Exported. " << "\n" << std::endl;
//      std::cout << "Matrix B Exported. " << "\n" << std::endl;


#else
  // Avoid compiler warnings
  libmesh_ignore(es);
#endif // LIBMESH_HAVE_SLEPC
}

void write_petsc_eigenvector_ascii(const EPS& eps,
                                   PetscInt eig_index,
                                   const libMesh::EquationSystems& equation_systems,
                                   const std::string& system_name,
                                   const std::string& var_name,
                                   const std::string& output_filename)
{
    // Get mesh and system info
    const libMesh::MeshBase& mesh = equation_systems.get_mesh();
    const libMesh::System& system = equation_systems.get_system(system_name);
    const libMesh::DofMap& dof_map = system.get_dof_map();
    const unsigned int var_num = system.variable_number(var_name);

    // Get the eigenvector from PETSc
    Vec xr;
    auto& libmesh_vec = dynamic_cast<libMesh::PetscVector<Number>&>(
		    *equation_systems.get_system("Eigensystem").solution
    );
    VecDuplicate(libmesh_vec.vec(), &xr);
    EPSGetEigenvector(eps, eig_index, xr, nullptr);

    // Open output file
    std::ofstream out(output_filename);
    out << std::scientific << std::setprecision(10);

    std::vector<libMesh::dof_id_type> dof_indices;
    for (const auto& node : mesh.local_node_ptr_range())
    {
        dof_map.dof_indices(node, dof_indices, var_num);
        if (dof_indices.empty())
            continue;

        PetscScalar psi_val = 0.0;
        PetscErrorCode ierr = VecGetValues(xr, 1, &dof_indices[0], &psi_val);
        if (ierr != 0) continue;

        out << (*node)(0) << " "  // x
            << (*node)(1) << " "  // y
            << (*node)(2) << " "  // z
            << PetscRealPart(psi_val) << "\n";
    }

    VecDestroy(&xr);
}

