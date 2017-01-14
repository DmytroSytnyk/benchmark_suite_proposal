#include "LBField.h"
#include <cstdlib>


    const int LBField::C[9][2] = { {-1,-1}, { 0,-1}, { 1,-1},
                                   {-1, 0}, { 0, 0}, { 1, 0},
                                   {-1, 1}, { 0, 1}, { 1, 1} };

    const double LBField::W[9] = { 1.0/36.0, 1.0/9.0, 1.0/36.0,
                                   1.0/9.0,  4.0/9.0, 1.0/9.0,
                                   1.0/36.0, 1.0/9.0, 1.0/36.0 };

    LBField::LBField(int nx, int ny):
    _nx(nx), _ny(ny),
    _pdf1(new double[(nx+2)*(ny+2)*9]), _pdf2(new double [(nx+2)*(ny+2)*9]),
    _pdfCollide(_pdf1), _pdfStream(_pdf2),
    _flag(new LBField::Flag[(nx+2)*(ny+2)]){
      if (_pdf1==NULL){std::cout << "ERROR LBField: _pdf1==NULL!" << std::endl; exit(EXIT_FAILURE);}
      if (_pdf2==NULL){std::cout << "ERROR LBField: _pdf2==NULL!" << std::endl; exit(EXIT_FAILURE);}
      if (_flag==NULL){std::cout << "ERROR LBField: _flag==NULL!" << std::endl; exit(EXIT_FAILURE);}
    }


    LBField::~LBField(){
      delete [] _pdf1; _pdf1=NULL;
      delete [] _pdf2; _pdf2=NULL;
      _pdfCollide=NULL; _pdfStream=NULL;
      delete [] _flag; _flag=NULL;
    }


    int LBField::getNx() const { return _nx;}
    int LBField::getNy() const { return _ny;}

    int LBField::getIndex(int x, int y) const { return (x+(_nx+2)*y); }

    double* LBField::getPdfCollide(){ return _pdfCollide;}
    double* LBField::getPdfStream() { return _pdfStream;}
    LBField::Flag*   LBField::getFlagField() { return _flag;}

    void LBField::swap(){
      double *tmp = _pdfCollide;
      _pdfCollide = _pdfStream;
      _pdfStream  = tmp;
    }

