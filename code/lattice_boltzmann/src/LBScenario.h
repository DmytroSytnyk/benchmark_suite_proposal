#ifndef _LBSCENARIO_H_
#define _LBSCENARIO_H_

#include "LBField.h"
#include <assert.h>
#include <string>
#include <cstdlib>
#include <cmath>

/** holds all parameters that are required to run a LB simulation.
 *  @author Philipp Neumann
 */
class LBParameters {
  public:
    int _nx; // number of inner LB nodes in x-direction
    int _ny; // number of inner LB nodes in y-direction
    int _timesteps; // number of LB time steps
    int _plotTimestep; // every interval of "plotTimestep": write a VTK
    double _tau;  // BGK relaxation time (in LB units)
    double _characteristicU; // characteristic velocity (in LB units)
    std::string _name; // name of simulated scenario
    double _positionX; // position of spherical obstacle in x-direction (scenario channel-obstacle)
    double _positionY; // position of spherical obstacle in y-direction (scenario channel-obstacle)
    double _radius;    // radius of spherical obstacle (scenario channel-obstacle)
};

/** implements difference flow scenarios.
 *  @author Philipp Neumann
 */
class LBScenario {
  public:
    LBScenario(){}
    virtual ~LBScenario(){}

    void initialiseField(LBField &field){
      const int nx = field.getNx();
      const int ny = field.getNy();
      double *pdfC = field.getPdfCollide();
      double *pdfS = field.getPdfStream();
      LBField::Flag   *flag = field.getFlagField();

      int index=0;
      for (int y = 0; y < ny+2; y++){
        for (int x = 0; x < nx+2; x++){
          initialiseCell(x,y,index,pdfC,pdfS,flag,nx,ny);
          index++;
        }
      }
    }

    /** initialises pdf values and flag field for a simulation. */
    virtual void initialiseCell(int x,int y, int index, double *pdfCollide, double *pdfStream, LBField::Flag* flag,const int nx, const int ny) = 0;
    /** returns true, if the cell x,y has a moving wall velocity and in this case sets ux, uy respectively. Returns false otherwise. */
    virtual bool getMovingWallVelocity(int x,int y, double &ux,double &uy,const int nx,const int ny) = 0;
};

// different initialisations ========================

/** initialisation without any physical boundary conditions anywhere. This only makes sense together with the LBNS implementation, where the
 *  first non-ghost LB nodes are overwritten each time step with boundary conditions from NS.
 *  @author Philipp Neumann
 */
class NoBoundary: public LBScenario {
  public:
    NoBoundary(){}
    virtual ~NoBoundary(){}

    virtual bool getMovingWallVelocity(int x, int y,double &ux,double& uy, const int nx, const int ny){ return false;}
    virtual void initialiseCell(int x, int y, int index, double *pdfCollide, double *pdfStream, LBField::Flag* flag, const int nx,const int ny){
      for (int i = 0; i < 9; i++){ pdfCollide[9*index+i] = LBField::W[i]; pdfStream[9*index+i] = LBField::W[i]; }
      flag[index] = LBField::FLUID;
    }
};

/** lid-driven cavity implementation
 *  @author Philipp Neumann
 */
class LBLidDrivenCavity: public LBScenario {
  public:
    LBLidDrivenCavity(double uxWall): _uxWall(uxWall){}
    virtual ~LBLidDrivenCavity(){}

    virtual bool getMovingWallVelocity(int x,int y, double &ux,double &uy,const int nx,const int ny){
      if (y>ny){ ux=_uxWall; uy=0.0; return true; }
      else     { return false; }
    }
    virtual void initialiseCell(int x,int y, int index,double *pdfCollide, double *pdfStream, LBField::Flag* flag,const int nx, const int ny){
      for (int i = 0; i < 9; i++){ pdfCollide[9*index+i] = LBField::W[i]; pdfStream[9*index+i] = LBField::W[i]; }
      flag[index] = LBField::FLUID;
      if (y > ny){ flag[index]=LBField::MOVINGWALL; }
      else {
        if ( (x<1) || (x>nx) || (y<1) ) { flag[index]=LBField::NOSLIP; }
      }

      // testing
      if ( (x>0.8*nx) && (x<0.9*nx) && (y>0.8*ny) && (y<0.9*ny)){flag[index]=LBField::NOSLIP;}
    }

  private:
    const double _uxWall; // velocity of lid in x-direction
};

/** couette flow implementation. The lower wall is being moved at constant velocity V. The upper wall is located at distance H.
 *  The left and right boundaries are periodic.
 *  The kinematic viscosity is given by nu.
 *  The analytic solution of the problem is given by:
 *  u(y,t)= V(1-y/H) - 2V/pi*sum_k=1^infty 1/k*sin(k*pi*y/H)*exp(-k^2 pi^2/H^2 * nu*t)
 *  @author Philipp Neumann
 */
class LBCouetteFlow: public LBScenario {
  public:
    LBCouetteFlow(double uxWall): _uxWall(uxWall){}
    virtual ~LBCouetteFlow(){}

    virtual bool getMovingWallVelocity(int x,int y, double &ux,double &uy,const int nx,const int ny){
      if (y<1){ ux=_uxWall; uy=0.0; return true; }
      else     { return false; }
    }
    virtual void initialiseCell(int x,int y, int index,double *pdfCollide, double *pdfStream, LBField::Flag* flag,const int nx, const int ny){
      for (int i = 0; i < 9; i++){ pdfCollide[9*index+i] = LBField::W[i]; pdfStream[9*index+i] = LBField::W[i]; }
      flag[index] = LBField::FLUID;
      if (y < 1){ flag[index]=LBField::MOVINGWALL; }
      else if (y > ny){ flag[index]=LBField::NOSLIP; }
      else {
        if ( (x<1) || (x>nx) ) { flag[index]=LBField::PERIODIC; }
      }
    }

  private:
    const double _uxWall; // velocity of lid in x-direction
};


/** taylor-green decaying vortices; we choose the form from wikipedia:
 *  - the velocity is assumed to be 1.0; characteristicU is the respective velocity in LB units.
 *  - domain size is 2*pi. This yields a meshsize dx=2.0*pi/nx.
 *  - the time step dt arises from the velocities and dx, see below.
 *  - using dx,dt, we rescale the initial solution to the LB velocities. The initial solution is:
 *    u(x,y)=U*sin(x)*cos(y), v(x,y)=-U*cos(x)*sin(y)
 *  - the initial pressure is:
 *    p(x,y)=-U*U/4*(cos(2x)+cos(2y))
 */
class LBTaylorGreen: public LBScenario {
  public:
    LBTaylorGreen(double characteristicU,double tau): _characteristicU(characteristicU), _tau(tau){}
    virtual ~LBTaylorGreen(){}

    virtual bool getMovingWallVelocity(int x,int y, double &ux,double &uy,const int nx,const int ny){return false;}
    virtual void initialiseCell(int x,int y, int index,double *pdfCollide, double *pdfStream, LBField::Flag* flag,const int nx, const int ny){
      const double posX= 2.0*3.141592653589793238*(((double)x)-0.5)/nx;
      const double posY= 2.0*3.141592653589793238*(((double)y)-0.5)/ny;

      const double ux = _characteristicU*cos(posX)*sin(posY);
      const double uy = -_characteristicU*sin(posX)*cos(posY);
      const double u2  = ux*ux+uy*uy;
      // stresses of form 0.5*(du_dy+dv_dx)
      const double sxx = _characteristicU*(-sin(posX))*sin(posY);
      const double syy = -_characteristicU*sin(posX)*(-sin(posY));
      const double sxy = 0.5*(_characteristicU*cos(posX)*cos(posY) - _characteristicU*cos(posX)*cos(posY));

      // shift by offset of default pressure
      const double p =1.0/3.0 -_characteristicU*_characteristicU/4.0*(cos(2*posX) + cos(2*posY));
      const double density=3.0*p;

      for (int i = 0; i < 9; i++){
        // compute equilibrium part from density and velocity
        const double cu=ux*LBField::C[i][0] + uy*LBField::C[i][1];
        const double feq = LBField::W[i]*density*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*u2);
        // Latt-based regularization to initialize non-equilibrium parts with stresses
        double fneq= LBField::W[i]*density*3.0*(-_tau)* (  (LBField::C[i][0]*LBField::C[i][0]-1.0/3.0)*sxx
                               + (LBField::C[i][1]*LBField::C[i][1]-1.0/3.0)*syy
                               +2.0       *LBField::C[i][0]*LBField::C[i][1]*sxy);
        pdfCollide[9*index+i] = feq+fneq;
        pdfStream[9*index+i]  = feq+fneq;
      }

      flag[index] = LBField::FLUID;
      if ( (y > ny) || (y <1) || (x>nx) || (x<1)){ flag[index]=LBField::PERIODIC; }
    }

  private:
    const double _characteristicU; // magnitude of vortex at beginning
    const double _tau; // relaxation time
};

class LBChannelWithObstacle: public LBScenario {
  public:
    LBChannelWithObstacle(double characteristicU, double viscosity, int nx, int ny, double posX, double posY, double radius):
    _characteristicU(characteristicU), _viscosity(viscosity), _Lx(nx*1.0), _Ly(ny*1.0), _posX(posX), _posY(posY), _radius(radius){
      assert (_posX-_radius>1.0); // obstacle at least one cell away from inlet
      assert (_Lx-(_posX+_radius) > 1.0); // obstacle at least one cell away from outlet
    }
    virtual ~LBChannelWithObstacle(){}

    virtual bool getMovingWallVelocity(int x,int y, double &ux,double &uy,const int nx,const int ny){
      // if this is a cell inside the left wall AND not part of the upper/lower noslip walls: return true and set parabolic vel. shape
      if ( (x < 1) && (y>0) && (y<ny+1)){
        const double posY = ((double)y)-0.5; // TODO check if this is correct; do we want to have the velocity of the neighbouring cell, or at the wall? Does it matter?
        uy=0.0;
        ux=4.0*_characteristicU*posY*(_Ly-posY)/(_Ly*_Ly);
        return true;
      } else { return false;}
    }
    virtual void initialiseCell(int x,int y, int index,double *pdfCollide, double *pdfStream, LBField::Flag* flag,const int nx, const int ny){
      for (int i = 0; i < 9; i++){pdfCollide[9*index+i]=LBField::W[i]; pdfStream[9*index+i]=LBField::W[i];}
      // lower/upper boundaries: noslip
      if ( (y < 1) || (y > ny) ){ flag[index]=LBField::NOSLIP;}
      // left: inlet
      else if (x<1){ flag[index]=LBField::MOVINGWALL;}
      // right: outlet
      else if (x>nx){flag[index]=LBField::PRESSURE;}
      // rest: either obstacle or fluid
      else {
        const double midPointX= ((double)x)-0.5;
        const double midPointY= ((double)y)-0.5;
        const double radius2  =_radius*_radius;
        if ( (midPointX-_posX)*(midPointX-_posX)+(midPointY-_posY)*(midPointY-_posY) < radius2 ){ flag[index]=LBField::NOSLIP;}
        else { flag[index]=LBField::FLUID;}
      }
    }

  private:
    const double _characteristicU;
    const double _viscosity;
    const double _Lx;
    const double _Ly;
    const double _posX;
    const double _posY;
    const double _radius;
};


class LBScenarioFactory {
  public:
  static LBScenario* getScenario(const LBParameters& parameters){
    if (parameters._name=="cavity"){ return new LBLidDrivenCavity(parameters._characteristicU);}
    else if (parameters._name=="couette"){return new LBCouetteFlow(parameters._characteristicU);}
    else if (parameters._name=="taylorgreen"){ return new LBTaylorGreen(parameters._characteristicU,parameters._tau);}
    else if (parameters._name=="channel-obstacle"){ return new LBChannelWithObstacle(
                                          parameters._characteristicU,(parameters._tau-0.5)/3.0,parameters._nx,parameters._ny,
                                          parameters._positionX,parameters._positionY,parameters._radius);}
    else if (parameters._name=="no-boundary"){ return new NoBoundary(); }
    else {std::cout << "ERROR LBScenarioFactory: Unknown scenario " << parameters._name << "!" << std::endl; exit(EXIT_FAILURE); return NULL; }
  }
};


#endif //
