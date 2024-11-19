// Revised version works with 2D PBC; confirmed using free rotor PES and G-matrix from 
// 1,4-dihydroxybenzene.  FEMvibn predicts spectrum with B value within 1 cm-1 of my
// estimate (21.7 vs 20.8 cm-1) and correct degeneracy.
// ALC 28-May-2024
// This version is working with 2D and 3D, accepting ndim and ElemType as arguments
// Not working with 1D -- throws memory access errors from PETSC, looks like something
//   being called that's not equipped for 1D matrices.
//
// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// <h1>Eigenproblems Example 2 - Solving a generalized Eigen Problem</h1>
// \author Steffen Petersen
// \date 2006
//
// This example shows how the previous EigenSolver example
// can be adapted to solve generalized eigenvalue problems.
//
// For solving eigen problems, libMesh interfaces
// SLEPc (www.grycap.upv.es/slepc/) which again is based on PETSc.
// Hence, this example will only work if the library is compiled
// with SLEPc support enabled.
//
// In this example some eigenvalues for a generalized symmetric
// eigenvalue problem A*x=lambda*B*x are computed, where the
// matrices A and B are assembled according to stiffness and
// mass matrix, respectively.

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
#include <iostream>
#include <string>
#include <unordered_map>
#include <getopt.h>
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
  // Initialize libMesh and the dependent libraries.
  LibMeshInit init (argc, argv);

  // This example is designed for the SLEPc eigen solver interface.
#ifndef LIBMESH_HAVE_SLEPC
  if (init.comm().rank() == 0)
    libMesh::err << "ERROR: This example requires libMesh to be\n"
                 << "compiled with SLEPc eigen solvers support!"
                 << std::endl;
#else

#ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
  // SLEPc currently gives us a nasty crash with Real==float
  libmesh_example_requires(false, "--disable-singleprecision");
#endif

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
//   d=ndim f=input_file e=nev p=npr q=npr1 t=elemTypeStr 
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
        std::cerr << "Usage: " << argv[0] << " [-f input_file] [-d ndim] [-e nev] [-p npr] [-q npr1] [t elemTypeStr] [-x nx] [-y ny] [-z nz] [-0 x0] [-1 x1] [-2 y0] [-3 y1] [-4 z0] [-5 z1] [-6 pbc_x] [-7 pbc_y] [-8 pbc_z]" << std::endl;
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
    std::cout << "elemTypeStr, ndim, nev, nx = " << elemTypeStr << " " << ndim << " " << nev << " " << nx << "\n";
  if (ndim == 2)
    std::cout << "elemTypeStr, ndim, nev, nx, ny = " << elemTypeStr << " " << ndim << " " << nev << " " << nx << " " << ny << "\n";
  if (ndim == 3)
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

        std::cout << "periodic boundary conditions flags:" << pbc_x << " " 
	          << pbc_y << " " << pbc_z << std::endl;
  // PERIODIC BOUNDARY CONDITIONS 
    DofMap & dof_map = eigen_system.get_dof_map();
    if ((ndim == 1) && (pbc_x == 1)){
        std::cout << "periodic boundary conditions for x" << std::endl;
	PeriodicBoundary xbdry(RealVectorValue(x1-x0));
        xbdry.myboundary = 0;
        xbdry.pairedboundary = 1;
        eigen_system.get_dof_map().add_periodic_boundary(xbdry);
    }

    if ((ndim == 2) && ((pbc_x == 1) || (pbc_y == 1))){
        std::cout << "periodic boundary conditions for " ;

	if ((pbc_x == 1) && (pbc_y == 1)){
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
          std::cout << "y" << std::endl;

          PeriodicBoundary xbdry(RealVectorValue(x1-x0, 0.0)); 
          xbdry.myboundary = 3;
          xbdry.pairedboundary = 1;
          eigen_system.get_dof_map().add_periodic_boundary(xbdry); 
	}
	if ((pbc_x == 1) && (pbc_y == 0)){
          std::cout << "x" << std::endl;

          PeriodicBoundary ybdry(RealVectorValue(0.0, y1-y0));
	  ybdry.myboundary = 0;
          ybdry.pairedboundary = 2;
          eigen_system.get_dof_map().add_periodic_boundary(ybdry); 
	}
    }

    if ((ndim == 3) && ((pbc_x == 1) || (pbc_y == 1) || (pbc_z == 1))){
        std::cout << "periodic boundary conditions for " ;

	if ((pbc_x == 1) && (pbc_y == 1) && (pbc_z == 1)){
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
          std::cout << "x" << std::endl;

          PeriodicBoundary xbdry(RealVectorValue(x1-x0, 0.0, 0.0)); 
          xbdry.myboundary = 4;
          xbdry.pairedboundary = 2;
          eigen_system.get_dof_map().add_periodic_boundary(xbdry); 
	}
	if ((pbc_x == 0) && (pbc_y == 1) && (pbc_z == 0)){
          std::cout << "y" << std::endl;

          PeriodicBoundary ybdry(RealVectorValue(0.0, y1-y0, 0.0));
	  ybdry.myboundary = 1;
          ybdry.pairedboundary = 3;
          eigen_system.get_dof_map().add_periodic_boundary(ybdry); 
	}
	if ((pbc_x == 0) && (pbc_y == 0) && (pbc_z == 1)){
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
  mesh.print_info();

        mesh_refinement.set_periodic_boundaries_ptr(dof_map.get_periodic_boundaries());


int n_nodes = mesh.n_nodes();
int n_elem = mesh.n_elem();
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

std::cout << "Writing number of dimensions to file: results/ndim" << std::endl;
  std::ofstream ndimfile;
  ndimfile.open("results/ndim");
  ndimfile << ndim ;
  ndimfile.close();

  // Declare the system variables.
  // Adds the variable "p" to "Eigensystem".   "p"
  // will be approximated using second-order approximation.
  eigen_system.add_variable("p", FIRST);

  // Give the system a pointer to the matrix assembly
  // function defined below.
  eigen_system.attach_assemble_function (assemble_mass);

  // Set necessary parameters used in EigenSystem::solve(),
  // i.e. the number of requested eigenpairs nev and the number
  // of basis vectors ncv used in the solution algorithm. Note that
  // ncv >= nev must hold and ncv >= 2*nev is recommended.
  equation_systems.parameters.set<unsigned int>("eigenpairs")    = nev;
  equation_systems.parameters.set<unsigned int>("basis vectors") = nev*3;

  // You may optionally change the default eigensolver used by SLEPc.
  // The Krylov-Schur method is mathematically equivalent to implicitly
  // restarted Arnoldi, the method of Arpack, so there is currently no
  // point in using SLEPc with Arpack.
  // ARNOLDI     = default in SLEPc 2.3.1 and earlier
  // KRYLOVSCHUR default in SLEPc 2.3.2 and later
  // eigen_system.get_eigen_solver().set_eigensolver_type(KRYLOVSCHUR);

// femvib-specific
      eigen_system.eigen_solver->set_eigensolver_type(KRYLOVSCHUR);
      eigen_system.eigen_solver->set_position_of_spectrum(SMALLEST_REAL);
// femvib-specific

  // Set the solver tolerance and the maximum number of iterations.
  equation_systems.parameters.set<Real>("linear solver tolerance") = pow(TOLERANCE, 5./3.);
  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 1000;

  // Set the type of the problem, here we deal with
  // a generalized Hermitian problem.
  eigen_system.set_eigenproblem_type(GHEP);

  // Set the eigenvalues to be computed. Note that not
  // all solvers in SLEPc support this capability.
//  eigen_system.get_eigen_solver().set_position_of_spectrum(2.3);
  ////  ALC commented out above line

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Prints information about the system to the screen.
  equation_systems.print_info();

#if SLEPC_VERSION_LESS_THAN(3,1,0)
  libmesh_error_msg("SLEPc 3.1 is required to call EigenSolver::set_initial_space()");
#else
  // Get the SLEPc solver object and set initial guess for one basis vector
  // this has to be done _after_ the EquationSystems object is initialized
  EigenSolver<Number> & slepc_eps = eigen_system.get_eigen_solver();
  NumericVector<Number> & initial_space = eigen_system.add_vector("initial_space");
  initial_space.add(1.0);
  slepc_eps.set_initial_space(initial_space);
#endif

  // Solve the system "Eigensystem".
  eigen_system.solve();

  // Get the number of converged eigen pairs.
  unsigned int nconv = eigen_system.get_n_converged();

  libMesh::out << "Number of converged eigenpairs: "
               << nconv
               << "\n"
               << std::endl;

  // Print the requested eigentates to ExodusII format files
  if (nconv != 0)
    {
#ifdef LIBMESH_HAVE_EXODUS_API
      for (int i=npr1; i<npr1+npr; i++)
      {
        std::ostringstream eigenstate_output_name;
        eigenstate_output_name << "results/eigenstate" << i << ".e";
        eigen_system.get_eigenpair(i);
        ExodusII_IO (mesh).write_equation_systems (eigenstate_output_name.str(), equation_systems);
      //      The gnuplot commands below may be useful for 1-D
//      GnuPlotIO plot(mesh, "Adaptivity Example 1", GnuPlotIO::GRID_ON);
//      plot.write_equation_systems ("gnuplot_script", equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
      }
      // next lines for extracting eigenvalues adapted from eigenproblem_ex3.C
        // write out all of the computed eigenvalues and plot the specified eigenvector
      std::ostringstream eigenvalue_output_name;
      eigenvalue_output_name << "results/eigenvalues.txt";
      std::ofstream evals_file(eigenvalue_output_name.str().c_str());
      for (unsigned int i=0; i<nconv; i++)
        {
          std::pair<Real,Real> eval = eigen_system.get_eigenpair(i);
          // The eigenvalues should be real!
          libmesh_assert_less (eval.second, TOLERANCE);
          evals_file << eval.first << std::endl;
        }
      evals_file.close();
    }
  else
    {
      libMesh::out << "WARNING: Solver did not converge!\n" << nconv << std::endl;
    }

#endif // LIBMESH_HAVE_SLEPC

  // All done.
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

std::cout << "Reading Kriging parameters from files results/optr results/opp" << std::endl;
  std::ifstream oppmfile;
  oppmfile.open("results/opp");
  oppmfile >> nptsx;
  oppmfile.close();
  std::ifstream optrfile;
  optrfile.open("results/optr");
  optrfile >> rval;
  optrfile.close();
 std::cout << "rval = " << rval << "\n" ;
 std::cout << "nptsx = " << nptsx << "\n" ;

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
  std::cout << "Writing PES data to file: results/pot.out " << std::endl;

// HERE IS WHERE WE WILL CHECK FOR DUPLICATES
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
      std::cout << "Finished writing PES data to file: " << "results/pot.out" << std::endl;
      std::cout << "Matrix A and B assembly completed. " << "\n" << std::endl;

      matrix_A.print_matlab("results/A.m");
      matrix_B.print_matlab("results/B.m");
      std::cout << "Matrix A Exported. " << "\n" << std::endl;
      std::cout << "Matrix B Exported. " << "\n" << std::endl;


#else
  // Avoid compiler warnings
  libmesh_ignore(es);
#endif // LIBMESH_HAVE_SLEPC
}
