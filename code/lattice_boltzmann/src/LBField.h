#ifndef _LBFIELD_H_
#define _LBFIELD_H_
#include <iostream>


/** simple data structure for D2Q9-based LB simulation.
 *  @author Philipp Neumann
 */
class LBField{
  public:
    enum Flag{FLUID=0,NOSLIP=1,MOVINGWALL=2,PERIODIC=3,PRESSURE=4};
    LBField(int nx, int ny);
    ~LBField();

    /** return domain size in x- and y-direction */
    int getNx() const;
    int getNy() const;

    /** returns the linearized index of a cell in the respective array */
    int getIndex(int x, int y) const;

    /** access to the data structures */
    double* getPdfCollide();
    double* getPdfStream();
    Flag*   getFlagField();

    /** const variant */
    const double * const getPdfCollide() const { return _pdfCollide;}

    /** swap access to pdf fields */
    void swap();

    static const int C[9][2]; // lattice velocities
    static const double W[9]; // lattice weights

  private:
    const int _nx; // size in x-direction
    const int _ny; // size in y-direction
    double *_pdf1; // pdf-field 1
    double *_pdf2; // pdf-field 2
    double *_pdfCollide; // (pre-)collide pdf field
    double *_pdfStream;  // pdf field that data are streamed to
    Flag*  _flag;
};

#endif // LBFIELD
