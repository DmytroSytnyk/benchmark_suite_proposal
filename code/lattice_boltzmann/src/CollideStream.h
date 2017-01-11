#ifndef _COLLIDESTREAM_H_
#define _COLLIDESTREAM_H_
#include "LBField.h"
#include "LBScenario.h"

/** implements a LB time step incl. boundary treatment.
 *  @author Philipp Neumann
 */
class CollideStream{
  public:
    CollideStream(LBScenario &scenario): _scenario(scenario){}
    virtual ~CollideStream(){}
    /** carries out collide stream and needs to call handleBoundary(...) for each cell. */
    virtual void collideStream(LBField &field) const = 0;

  protected:
    /** applies boundary treatment to the cells. */
    void handleBoundary(
      const int x,const int y, const int cellIndex, const int nbX, const int nbY, const int nbCell, const int pdfIndex, LBField& field,
      const double postCollisionPdf, const double cellDensity
    ) const;

  private:
    LBScenario &_scenario;
};


/** implements the BGK dynamics. */
class BGK: public CollideStream {
  public:
    BGK(LBScenario &scenario, double tau): CollideStream(scenario), _tau(tau), _omega(1.0/tau){}
    virtual ~BGK(){}

    void collideStream(LBField& field) const;

    /** carries out normalisation of pdfs in collide field right before collide step, following the procedure in
     *  "Lattice Boltzmann method with regularized pre-collision distribution functions" by J. Latt+B. Chopard.
     */
    void normalisePdfsLatt(LBField &field) const ;
    /** carries out normalisation of pdfs based on Neumann approach, cf. "Hybrid Multiscale Simulation Approaches for Micro- and Nanoflows" by P. Neumann */
    void normalisePdfsNeumann(LBField &field) const;
private:
    const double _tau;  // relaxation time
    const double _omega;// relaxation frequency (1.0/_tau)
};

#endif // COLLIDESTREAM
