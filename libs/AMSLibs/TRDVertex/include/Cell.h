#ifndef CELL_H
#define CELL_H

#include "TRD2DHit.h"
#include <vector>
#include <TChain.h>
#include <iostream>
#include <TMath.h>
#include <TString.h>

class TRD2DHit;

class Cell{

public:
    Cell();
    Cell(TRD2DHit* _p1,TRD2DHit*_p2);
    ~Cell();

    bool IsConnected(Cell* cell2);
    bool IsNeighbour(Cell* cell2);
    void Print();
    void PrintLineRepresentation();
    void Draw();
    bool operator() (Cell i,Cell j) { return (i.state>j.state);}

    void AddConnectedCell();
    void ClearConnectedCell();

    void Evolve();

    double GetTransformFunction(Cell* cell2,double alpha,double gamma);

    bool IsSameState(Cell *cell2);
    bool IsMatched(Cell *cell2);

    void InValidate();

    void InValidateOtherCells();

    bool IsValid(){return valid && state>0;}
    void SortConnectedCells();

    int Merge(Cell *cell2);
    double CosTheta(Cell *cell2);

    //    double GetAngle(Cell *cell2);
    //    void GetAngleRange(Cell *cell2, double *range);
    //    double GetAbsoluteAngle(Cell *cell2);
    bool Intersect(Cell *c2, float &i_x, float &i_y);





    bool valid;

    bool HasMatchedCells;
    TRD2DHit* p1;
    TRD2DHit* p2;
    double theta;
    double theta_range[2];
    double state;
    double evolve_next_generation;
    int assigned;

  //  std::vector<TRD2DHit*> v_hits;

    std::vector<Cell*> v_connected_upper;
    std::vector<Cell*> v_connected_lower;
    std::vector<Cell*> v_connected_samestart;
    std::vector<Cell*> v_connected_sameend;

    double connected_upper_minangle;
    double connected_lower_minangle;

    double length;

};


class compCellByCosTheta{
public:
    compCellByCosTheta(Cell* m):mother(m){}

    bool operator()( Cell *a,  Cell *b)const{
        return mother->CosTheta(a) > mother->CosTheta(b);
    }

    Cell* mother;

};


#endif // CELL_H
