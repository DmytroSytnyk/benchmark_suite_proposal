#ifndef _MAIN_CPP_
#define _MAIN_CPP_

#include "LBField.h"
#include "LBScenario.h"
#include "CollideStream.h"
#include "VTK.h"
#include <cstdlib>

/**
 * @author Philipp Neumann
 */
int main(int argc, char *argv[]){
  LBParameters parameters;
  parameters._nx = 32;         // number of grid cells x-direction
  parameters._ny = 32;         // number of grid cells y-direction
  parameters._timesteps=10000; // number of time steps
  parameters._plotTimestep=10; // plot every 10 time steps

  parameters._tau=0.65;             // relaxation time tau= kin.viscosity*3 + 0.5
  parameters._characteristicU=0.05; // characteristic velocity
  parameters._name="couette";   // for other scenarios: change this parameter, e.g. to "couette"
  parameters._positionX=parameters._nx/11.0;
  parameters._positionY=parameters._positionX;
  parameters._radius   =0.25*parameters._positionX;


  LBField field(parameters._nx,parameters._ny);
  LBScenario *scenario = LBScenarioFactory::getScenario(parameters);
  if (scenario==NULL){std::cout << "ERROR: scenario==NULL!" << std::endl; exit(EXIT_FAILURE);}
  scenario->initialiseField(field);
  BGK collideStream(*scenario,parameters._tau);

  for (int t = 0; t < parameters._timesteps; t++){
    // just for current test: only start plotting after 125 000 time steps (corresponding to simulation time t_NS=5.0)
    //if (t>=125000){
    VTK::plot(field,t,parameters._name,parameters._plotTimestep,parameters._tau);
    //}
    // test: use Latt normal. procedure to stabilise flow
    //collideStream.normalisePdfsNeumann(field);
    collideStream.collideStream(field);
    std::cout << "Finish timestep " << t << "..." << std::endl;
  }

  delete scenario; scenario = NULL;
  return 0;
}

#endif 
