#include "epot_bicgstabsolver.hpp"
#include "particledatabase.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "particlediagplotter.hpp"
#include "geomplotter.hpp"
#include "config.h"
#include "trajectorydiagnostics.hpp"
#include "fstream"
#ifdef GTK3
#include "gtkplotter.hpp"
#endif


const double Te = 7.0;
const double Up = 7.0;
 double l = 8e-3;
 double r_pul = 4e-3;
 double r_ext = 6e-3;

bool solid1( double x, double y, double z )
{
    return( x <= 9e-3 && x >= 2e-3 && y  >= r_pul  + (1.0*x-3e-3)&& y >= r_pul);
}


bool solid2( double x, double y, double z )
{
    return( x >= 9e-3 + l && y >= r_ext && y >= -1*(x- 9e-3 - l - 3e-3 ) + r_ext );
}


void simu( int *argc, char ***argv )
{
    Geometry geom( MODE_CYL, Int3D(500,241,1), Vec3D(0,0,0), 0.00005 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0 ) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, -40.0e3) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0)  );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, -40.0e3) );
    geom.build_mesh();

    EpotBiCGSTABSolver solver( geom );
    InitialPlasma initp( AXIS_X, 1e-3 );
    solver.set_initial_plasma( 7.0, &initp );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;
    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE,
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBaseCyl pdb( geom );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );

    for( size_t i = 0; i < 15; i++ ) {

	if( i == 1 ) {
	    double rhoe = pdb.get_rhosum();
	    solver.set_pexp_plasma( -rhoe, Te, Up );
	}

	solver.solve( epot, scharge );
	efield.recalculate();

	pdb.clear();
	pdb.add_2d_beam_with_energy( 15000, 1000.0, 1.0, 1.0, 
				     7.0, 0.0, 0.5, 
				     0.0, 0.0, 
				     0.0, 0.012 );
	pdb.iterate_trajectories( scharge, efield, bfield );
   }


TrajectoryDiagnosticData tdata;
std::vector<trajectory_diagnostic_e> diagnostics;
diagnostics.push_back( DIAG_Y );
diagnostics.push_back( DIAG_YP );
diagnostics.push_back(DIAG_CURR);
pdb.trajectories_at_plane( tdata, AXIS_X, geom.max(0)-geom.h(), diagnostics );
Emittance emit( tdata(0).data(), tdata(1).data());
double I = 0.0;
for(size_t i = 0; i < tdata.traj_size();i++){//sum current just behind the extractor
            I += tdata(2)[i];
    }
    

// Output
std::ofstream dout( "data/emit.txt", std::ios::app );
     dout << r_pul << " "
          << r_ext << " "
          << l << " "
          << emit.alpha() << " "
          << emit.beta() << " "
          << emit.epsilon() << " "
          << I << "/n";
     dout.close();

#ifdef GTK3
    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_scharge( &scharge );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
#endif
}


int main( int argc, char **argv )
{
    try {
        ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
        simu( &argc, &argv );
    } catch ( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}
