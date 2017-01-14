#include "CollideStream.h"

#include <cstdlib>
#include <assert.h>

void CollideStream::handleBoundary(
const int x,const int y, const int cellIndex, const int nbX, const int nbY, const int nbCell, const int pdfIndex,
LBField& field, const double postCollisionPdf, const double cellDensity) const{
   const LBField::Flag  * const flag = field.getFlagField();
   const int nx = field.getNx();
   const int ny = field.getNy();
   double *pdfS = field.getPdfStream();

  if (flag[nbCell]==LBField::NOSLIP){
    //std::cout << "Cell " << nbX << "," << nbY << " is noslip" << std::endl;
    pdfS[9*cellIndex+(8-pdfIndex)] = postCollisionPdf;
  }
  else if (flag[nbCell]==LBField::MOVINGWALL){
    //std::cout << "Cell " << nbX << "," << nbY << " is movingwall" << std::endl;
    double ux,uy;
    bool movingWall=_scenario.getMovingWallVelocity(nbX,nbY,ux,uy,nx,ny);
    if (!movingWall){std::cout << "ERROR CollideStream::handleBoundary: getMovingWallVelocity and FlagField disagree :-(!" << std::endl; exit(EXIT_FAILURE);}
    pdfS[9*cellIndex+(8-pdfIndex)] = postCollisionPdf - 6.0*LBField::W[pdfIndex]*cellDensity*(LBField::C[pdfIndex][0]*ux+LBField::C[pdfIndex][1]*uy);
  } else if (flag[nbCell]==LBField::PRESSURE){
    // use pressure condition (assume fixed pressure=unit density), e.g. see PhD thesis by Nils Thuerey, Eq. (4.5) 
    double *pdfCell = &(field.getPdfCollide()[9*cellIndex]); // access to pre-collide field
    double vel[2] ={0.0,0.0};                              // velocity in this cell
    for (int i = 0; i < 9; i++){ vel[0] += pdfCell[i]*LBField::C[i][0]; vel[1] +=pdfCell[i]*LBField::C[i][1]; }
    vel[0] = vel[0]/cellDensity; vel[1] = vel[1]/cellDensity;

    const double cu = vel[0]*LBField::C[pdfIndex][0]+vel[1]*LBField::C[pdfIndex][1];
    const double feq= LBField::W[pdfIndex]*1.0*(1.0+3.0*cu+4.5*cu*cu-1.5*(vel[0]*vel[0]+vel[1]*vel[1]));
    const double feqInv=feq - 6.0*1.0*LBField::W[pdfIndex]*cu;
    pdfS[9*cellIndex+(8-pdfIndex)] = feq+feqInv - pdfCell[pdfIndex];
  } else if (flag[nbCell]==LBField::PERIODIC){
    int targetX=nbX;
    int targetY=nbY;
    if (nbX==0){targetX=nx;} else if (nbX==nx+1){targetX=1;}
    if (nbY==0){targetY=ny;} else if (nbY==ny+1){targetY=1;}
    pdfS[9*field.getIndex(targetX,targetY)+pdfIndex] = postCollisionPdf;
  }
}


void BGK::collideStream(LBField& field) const {
  // debugging: carry out normalisation of pdfs to remove higher-order effects
  // normalisePdfsLatt(field);
 
  const int nx = field.getNx();
  const int ny = field.getNy();
  const double * const pdfC = field.getPdfCollide();
  double *pdfS = field.getPdfStream();
  const LBField::Flag   * const flag = field.getFlagField();

  int index = (nx+2)+1;
  int pdfIndex=9*index;

  for (int y=1; y < ny+1; y++){
    for (int x=1; x < nx+1; x++){
      // compute density+velocity
      double density=0.0;
      double velocity[2]={0.0,0.0};
      velocity[0] = pdfC[pdfIndex+2] + pdfC[pdfIndex+5]+pdfC[pdfIndex+8];
      density     = pdfC[pdfIndex]   + pdfC[pdfIndex+3]+pdfC[pdfIndex+6];
      velocity[1] = -pdfC[pdfIndex]  - pdfC[pdfIndex+1]-pdfC[pdfIndex+2]
                    +pdfC[pdfIndex+6]+ pdfC[pdfIndex+7]+pdfC[pdfIndex+8];
      velocity[0]-= density;
      density    +=  (pdfC[pdfIndex+1]+pdfC[pdfIndex+2]
                     +pdfC[pdfIndex+4]+pdfC[pdfIndex+5]
                     +pdfC[pdfIndex+7]+pdfC[pdfIndex+8]);
      velocity[0] = velocity[0]/density;
      velocity[1] = velocity[1]/density;
      const double u2= 1.0 - 1.5*(velocity[0]*velocity[0]+velocity[1]*velocity[1]);
      //std::cout << "BGK::collideStream: " << x << "," << y << ": " << density-1.0 << std::endl; assert(fabs(density-1.0)<1.0e-12); //debugging

      // compute equilibrium distribution, do BGK collide step+streaming =======
      const double w0 = LBField::W[0]*density;
      const double w1 = LBField::W[1]*density;
      // distribution 0,8
      double cu = LBField::C[0][0]*velocity[0]+LBField::C[0][1]*velocity[1];
      double feq= w0*(u2+3.0*cu+4.5*cu*cu);
      double postPdf = pdfC[pdfIndex] - _omega*(pdfC[pdfIndex]-feq);
      int nbX   = x+LBField::C[0][0]; int nbY = y+LBField::C[0][1];
      int nbCell= field.getIndex(nbX,nbY);
      if (flag[nbCell]==LBField::FLUID){
        // streaming step
        pdfS[9*nbCell] = postPdf;
      } else {
        handleBoundary(x,y,index,nbX,nbY,nbCell,0,field,postPdf,density);
      }
      feq     = feq - 6.0*w0*cu;
      postPdf = pdfC[pdfIndex+8] - _omega*(pdfC[pdfIndex+8]-feq);
      nbX   = x+LBField::C[8][0]; nbY = y+LBField::C[8][1];
      nbCell= field.getIndex(nbX,nbY);
      if (flag[nbCell]==LBField::FLUID){
        // streaming step
        pdfS[9*nbCell+8] = postPdf;
      } else {
        handleBoundary(x,y,index,nbX,nbY,nbCell,8,field,postPdf,density);
      }
      // distribution 1,7
      cu = LBField::C[1][1]*velocity[1];
      feq= w1*(u2+3.0*cu+4.5*cu*cu);
      postPdf = pdfC[pdfIndex+1] - _omega*(pdfC[pdfIndex+1]-feq);
      nbX   = x; nbY = y+LBField::C[1][1];
      nbCell= field.getIndex(nbX,nbY);
      if (flag[nbCell]==LBField::FLUID){
        // streaming step
        pdfS[9*nbCell+1] = postPdf;
      } else {
        handleBoundary(x,y,index,nbX,nbY,nbCell,1,field,postPdf,density);
      }
      feq     = feq - 6.0*w1*cu;
      postPdf = pdfC[pdfIndex+7] - _omega*(pdfC[pdfIndex+7]-feq);
      nbY = y+LBField::C[7][1];
      nbCell= field.getIndex(nbX,nbY);
      if (flag[nbCell]==LBField::FLUID){
        // streaming step
        pdfS[9*nbCell+7] = postPdf;
      } else {
        handleBoundary(x,y,index,nbX,nbY,nbCell,7,field,postPdf,density);
      }
      // distribution 2,6
      cu = LBField::C[2][0]*velocity[0]+LBField::C[2][1]*velocity[1];
      feq= w0*(u2+3.0*cu+4.5*cu*cu);
      postPdf = pdfC[pdfIndex+2] - _omega*(pdfC[pdfIndex+2]-feq);
      nbX   = x +LBField::C[2][0]; nbY = y+LBField::C[2][1];
      nbCell= field.getIndex(nbX,nbY);
      if (flag[nbCell]==LBField::FLUID){
        // streaming step
        pdfS[9*nbCell+2] = postPdf;
      } else { 
        handleBoundary(x,y,index,nbX,nbY,nbCell,2,field,postPdf,density);
      }
      feq     = feq - 6.0*w0*cu;
      postPdf = pdfC[pdfIndex+6] - _omega*(pdfC[pdfIndex+6]-feq);
      nbX   = x+LBField::C[6][0]; nbY = y+LBField::C[6][1];
      nbCell= field.getIndex(nbX,nbY);
      if (flag[nbCell]==LBField::FLUID){
        // streaming step
        pdfS[9*nbCell+6] = postPdf;
      } else {
        handleBoundary(x,y,index,nbX,nbY,nbCell,6,field,postPdf,density);
      }
      // distribution 3,5
      cu = LBField::C[3][0]*velocity[0];
      feq= w1*(u2+3.0*cu+4.5*cu*cu);
      postPdf = pdfC[pdfIndex+3] - _omega*(pdfC[pdfIndex+3]-feq);
      nbX   = x +LBField::C[3][0]; nbY = y;
      nbCell= field.getIndex(nbX,nbY);
      if (flag[nbCell]==LBField::FLUID){
        // streaming step
        pdfS[9*nbCell+3] = postPdf;
      } else { 
        handleBoundary(x,y,index,nbX,nbY,nbCell,3,field,postPdf,density);
      }
      feq     = feq - 6.0*w1*cu;
      postPdf = pdfC[pdfIndex+5] - _omega*(pdfC[pdfIndex+5]-feq);
      nbX   = x+LBField::C[5][0];
      nbCell= field.getIndex(nbX,nbY);
      if (flag[nbCell]==LBField::FLUID){
        // streaming step
        pdfS[9*nbCell+5] = postPdf;
      } else {
        handleBoundary(x,y,index,nbX,nbY,nbCell,5,field,postPdf,density);
      }
      // distribution 4
      feq = LBField::W[4]*density*u2;
      pdfS[pdfIndex+4]=pdfC[pdfIndex+4] - _omega*(pdfC[pdfIndex+4]-feq);

      // ================================================

      // update index (in x-line)
      index++;
      pdfIndex+=9;  
    }
    // update index (jumping into next x-line and skipping 2 ghost layer cells)
    index+=2;
    pdfIndex+=18;
  }

  field.swap();
}


void BGK::normalisePdfsLatt(LBField &field) const {
  const int nx = field.getNx();
  const int ny = field.getNy();
  double *pdfC = field.getPdfCollide();
  const LBField::Flag   * const flag = field.getFlagField();
  // for test: save lower left corner node and check on equality in LBNS simulation
  //double lowerLeft[9];
  //for (int i = 0; i < 9; i++){ lowerLeft[i] = pdfC[9*field.getIndex(1,1)+i]; }


  // loop over all inner cells x,y
  for (int y = 1; y < ny+1; y++){
    for (int x = 1; x < nx+1; x++){
      // only carry out normalisation on fluid nodes
      if (flag[field.getIndex(x,y)]==LBField::FLUID){
        double density = 0.0;
        double ux      = 0.0;
        double uy      = 0.0;
        double s_xx    = 0.0;
        double s_yy    = 0.0;
        double s_xy    = 0.0;
        double feq[9]  = {0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0};
        const int offset= 9*field.getIndex(x,y);
        // compute density+ velocity
        for (int i = 0; i < 9; i++){
          density += pdfC[offset+i];
          ux      += pdfC[offset+i]*LBField::C[i][0];
          uy      += pdfC[offset+i]*LBField::C[i][1];
        }
        ux = ux/density; uy = uy/density;
        // compute feq and stresses
        for (int i = 0; i < 9; i++){
          const double cu = ux*LBField::C[i][0] + uy*LBField::C[i][1];
          feq[i] = LBField::W[i]*density*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*(ux*ux+uy*uy));
          s_xx += (pdfC[offset+i]-feq[i])*LBField::C[i][0]*LBField::C[i][0];
          s_yy += (pdfC[offset+i]-feq[i])*LBField::C[i][1]*LBField::C[i][1];
          s_xy += (pdfC[offset+i]-feq[i])*LBField::C[i][0]*LBField::C[i][1];
        }
        // compute fneq and set distribution
        for (int i = 0; i < 9; i++){
          const double fneq = (  (LBField::C[i][0]*LBField::C[i][0]-1.0/3.0)*s_xx
                               + (LBField::C[i][1]*LBField::C[i][1]-1.0/3.0)*s_yy
                               +2.0       *LBField::C[i][0]*LBField::C[i][1]*s_xy)*LBField::W[i]*4.5;
          pdfC[offset+i] = feq[i]+fneq;
        }
      }
    }
  }

  // check pdfs on lower left corner
  /*for (int i = 0; i < 9; i++){
    assert(fabs(lowerLeft[i]-pdfC[9*field.getIndex(1,1)+i])<1.0e-14);
  }*/
}



void BGK::normalisePdfsNeumann(LBField &field) const {
  const int nx = field.getNx();
  const int ny = field.getNy();
  double *pdfC = field.getPdfCollide();
  const LBField::Flag   * const flag = field.getFlagField();

  const double optMatrix[9][3] = { { 1.0/36.0, 1.0/36.0, 0.25},
                                   {-1.0/18.0, 4.0/9.0 , 0.0 },
                                   { 1.0/36.0, 1.0/36.0,-0.25},
                                   { 4.0/9.0 ,-1.0/18.0, 0.0 },
                                   {-8.0/9.0 ,-8.0/9.0 , 0.0 },
                                   { 4.0/9.0 ,-1.0/18.0, 0.0 },
                                   { 1.0/36.0, 1.0/36.0,-0.25},
                                   {-1.0/18.0, 4.0/9.0 , 0.0 },
                                   { 1.0/36.0, 1.0/36.0, 0.25} };



  // loop over all inner cells x,y
  for (int y = 1; y < ny+1; y++){
    for (int x = 1; x < nx+1; x++){
      // only carry out normalisation on fluid nodes
      if (flag[field.getIndex(x,y)]==LBField::FLUID){
        double density = 0.0;
        double ux      = 0.0;
        double uy      = 0.0;
        double s_xx    = 0.0;
        double s_yy    = 0.0;
        double s_xy    = 0.0;
        double feq[9]  = {0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0};
        const int offset= 9*field.getIndex(x,y);
        // compute density+ velocity
        for (int i = 0; i < 9; i++){
          density += pdfC[offset+i];
          ux      += pdfC[offset+i]*LBField::C[i][0];
          uy      += pdfC[offset+i]*LBField::C[i][1];
        }
        ux = ux/density; uy = uy/density;
        // compute feq and stresses
        for (int i = 0; i < 9; i++){
          const double cu = ux*LBField::C[i][0] + uy*LBField::C[i][1];
          feq[i] = LBField::W[i]*density*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*(ux*ux+uy*uy));
          s_xx += (pdfC[offset+i]-feq[i])*LBField::C[i][0]*LBField::C[i][0];
          s_yy += (pdfC[offset+i]-feq[i])*LBField::C[i][1]*LBField::C[i][1];
          s_xy += (pdfC[offset+i]-feq[i])*LBField::C[i][0]*LBField::C[i][1];
        }
        // compute fneq and set distribution
        for (int i = 0; i < 9; i++){
          const double fneq = (optMatrix[i][0]*s_xx + optMatrix[i][1]*s_yy + optMatrix[i][2]*s_xy)/(_omega*3.0);
          pdfC[offset+i] = feq[i]+fneq;
        }
      }
    }
  }

  // check pdfs on lower left corner
  /*for (int i = 0; i < 9; i++){
    assert(fabs(lowerLeft[i]-pdfC[9*field.getIndex(1,1)+i])<1.0e-14);
  }*/
}
