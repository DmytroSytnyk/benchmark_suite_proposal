#ifndef _VTK_H_
#define _VTK_H_
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "LBField.h"

class VTK{
  public:
  /** plotting for BGK simulations: plots density, velocity, and stresses.
   *  If offsetX,offsetY are specified, the domain is shifted from the origin (0,0) for the lower left corner of the computational domain to offsetX,offsetY.
   *  If hx,dt are specified, a mesh size hx and time step dt are assumed. The velocity and stresses are scaled in this case respectively, following the formula
   *  velocity_VTK = velocity_LB * hx/dt, stress_VTK = stress_LB/dt
   *  field - LB field with pdfs
   *  t - time step to be written
   *  filestem - name of vtk file (will be extended by "LB_" in front, and "t.vtk" with t the time step)
   *  plotTimestep - number of LB time steps between two VTK plots
   *  tau - BGK relaxation time
   *
   *  @author Philipp Neumann
   */
  static void plot(LBField &field, int t, std::string filestem,int plotTimestep,double tau,double offsetX=0.0,double offsetY=0.0,double hx=1.0,double dt=1.0,bool plotGhostCell=true)  {
    // only plot output if this is the correct timestep
    if (plotTimestep==-1){ return;}
    if (t%plotTimestep!=0){return;}

    std::stringstream ss; ss << "LB_" << filestem << t << ".vtk";
    std::ofstream file(ss.str().c_str());
    if (!file.is_open()){std::cout << "ERROR VTK::plot(): Could not open file " << ss.str() << "!" << std::endl; exit(EXIT_FAILURE);}
    std::stringstream flag, density, velocity,s_xx_stream,s_xy_stream,s_yy_stream,pressure;
    LBField::Flag *flagField = field.getFlagField();
    const double* const pdfField  = field.getPdfCollide();

    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "Super-Epsi's simple LB simulation" << std::endl;
    file << "ASCII" << std::endl << std::endl;
    file << "DATASET STRUCTURED_GRID" << std::endl;
    file << "DIMENSIONS " << field.getNx()+1+2*plotGhostCell << " " << field.getNy()+1+2*plotGhostCell  << " " << 1 << std::endl;
    file << "POINTS " << (field.getNx()+1+2*plotGhostCell)*(field.getNy()+1+2*plotGhostCell) << " float" << std::endl;

    flag << "CELL_DATA " << (field.getNx()+2*plotGhostCell)*(field.getNy()+2*plotGhostCell) << std::endl;
    flag << "SCALARS flag float 1" << std::endl;
    flag << "LOOKUP_TABLE default" << std::endl;

    density  << std::setprecision(16);
    density  << "SCALARS density float 1 " << std::endl;
    density  << "LOOKUP_TABLE default" << std::endl;
    pressure << std::setprecision(16);
    pressure << "SCALARS pressure float 1 " << std::endl;
    pressure << "LOOKUP_TABLE default" << std::endl;

    s_xx_stream << std::setprecision(16); s_yy_stream << std::setprecision(16); s_xy_stream << std::setprecision(16);
    velocity << std::setprecision(16);
    s_xx_stream << "SCALARS s_xx float 1 " << std::endl; s_xx_stream << "LOOKUP_TABLE default" << std::endl;
    s_yy_stream << "SCALARS s_yy float 1 " << std::endl; s_yy_stream << "LOOKUP_TABLE default" << std::endl;
    s_xy_stream << "SCALARS s_xy float 1 " << std::endl; s_xy_stream << "LOOKUP_TABLE default" << std::endl;
    velocity << "VECTORS velocity float" << std::endl;

    // loop over domain (incl. boundary) and write point coordinates
    for (int y = -1*plotGhostCell; y < field.getNy()+1+plotGhostCell; y++){ for (int x = -1*plotGhostCell; x < field.getNx()+1+plotGhostCell; x++){
      file << x*hx+offsetX << " " <<  y*hx+offsetY  << " 0.0" << std::endl;
    }}

    // loop over domain (incl. boundary)
    for (int y = 1-plotGhostCell; y < field.getNy()+1+plotGhostCell; y++){ for (int x = 1-plotGhostCell; x < field.getNx()+1+plotGhostCell; x++){
      int index=field.getIndex(x,y);
      double dens = 1.0;
      double vel[2] = {0.0,0.0};
      double s_xx = 0.0; double s_xy = 0.0; double s_yy=0.0;
      // compute density and velocity for fluid cells only and set unit density/zero velocity otherwise
      if (flagField[index]==LBField::FLUID){
        dens=0.0;
        for (int i = 0; i < 9; i++){
          dens += pdfField[9*index+i];
          for (int d = 0; d < 2; d++){ vel[d] += pdfField[9*index+i]*LBField::C[i][d]; }
        }
        for (int d = 0; d < 2; d++){ vel[d] /= dens; }
        // compute stresses (s_xy = drho u_y/dx +  drho u_x/dy etc)
        for (int i = 0; i < 9; i++){
          const double cu= LBField::C[i][0]*vel[0]+LBField::C[i][1]*vel[1];
          const double feq = LBField::W[i]*dens*(1.0+3.0*cu+4.5*cu*cu-1.5*(vel[0]*vel[0]+vel[1]*vel[1]));
          s_xx += (pdfField[9*index+i]-feq)*LBField::C[i][0]*LBField::C[i][0];
          s_yy += (pdfField[9*index+i]-feq)*LBField::C[i][1]*LBField::C[i][1];
          s_xy += (pdfField[9*index+i]-feq)*LBField::C[i][0]*LBField::C[i][1];
        }
        s_xx = -s_xx/tau*3.0;
        s_yy = -s_yy/tau*3.0;
        s_xy = -s_xy/tau*3.0;
      }
      // write information to streams
      flag << flagField[index] << std::endl;
      density << dens << std::endl;
      pressure << (dens-1.0)/3.0*hx*hx/(dt*dt) << std::endl;
      s_xx_stream << s_xx/dt << std::endl; s_yy_stream << s_yy/dt << std::endl; s_xy_stream << s_xy/dt << std::endl;
      velocity << vel[0]*hx/dt << " " << vel[1]*hx/dt << " 0.0" << std::endl;
      index++;
    }}

    file << std::endl;
    file << flag.str() << std::endl << std::endl;
    flag.str("");
    file << density.str() << std::endl;
    density.str("");
    file << pressure.str() << std::endl;
    pressure.str("");
    file << s_xx_stream.str() << std::endl; s_xx_stream.str("");
    file << s_yy_stream.str() << std::endl; s_yy_stream.str("");
    file << s_xy_stream.str() << std::endl; s_xy_stream.str("");
    file << velocity.str() << std::endl;
    velocity.str("");
    file.close();
  }

};

#endif
