#ifndef _TRDVERTEX_
#define _TRDVERTEX_
#include "TF1.h"
#include "TF2.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "TMinuit.h"
#include "TGraph.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TLine.h"
#include <TMath.h>
#include <TString.h>
#include "Segment.h"
#include "Cell.h"
#include "Vertex2D.h"
#include "Vertex3D.h"
#include "TRD2DHit.h"
#include "HitState.h"
#include "TrdKHit.h"




bool comp2DVertexByNTracks(Vertex2D *a, Vertex2D *b);

bool comp3DVertexByNTracks(Vertex3D *a, Vertex3D *b);

bool comp2DVertexByNHits(Vertex2D *a, Vertex2D *b);

bool LessThanThreeHits(Line* line);

bool compByState(Cell* a, Cell* b);

bool compByLength(const Cell *a, const Cell *b);


class TRDVertex{

public:

    TRDVertex(){
        Init();
    }

    ~TRDVertex(){
        delete Collection_Tracks_xz;
        delete Collection_Tracks_yz;
        delete v_vertex;
        delete v_v2d_xz;
        delete v_v2d_yz;
    }

    void SetHit(std::vector<SimpleHitState> *MySimpleHits_TRD);
    void SetRandomHits();

    TRD2DHitCollection Hits;


    static bool DEBUG;

    void Init();

    int Reconstruction();
    int Reconstruction(AMSEventR* evt);
    void TrackFinding();

    void VertexFinding();

    void VertexFinding_New();

    std::vector<Cell*> cells_x;
    std::vector<Cell*> cells_y;

    void CleanUp(std::vector<Cell*>& cells);

    void Evolve(std::vector<Cell*>& cells, int iside);

    void SelectTracks(std::vector<Cell*>& cells, int iside);

    void SelectTracks_New(std::vector<Cell*>& cells, int iside);
    void MakeCell_New(std::vector<Cell*> &cells, TString Projection);

    void MakeCell_New_Improved(std::vector<Cell*> &cells, TString Projection);

    void Evolve_NN(std::vector<Cell*>& cells,int iside);

    void DrawVertex(int iside, Vertex3D *v3d, int index);
    void DrawVertex(int iside, Vertex2D *v2d, int index);
    void Build2DVertex(int iside);
    void DrawCells(std::vector<Cell*> &cells);


    void Build2DVertex_ZVTOP(int iside);
    void ClearMemory();
    void MyLabel(const char *drawstring, Double_t x, Double_t y, float size, Color_t color = kBlack);
    void Build2DVertex_Adaptive(int iside);

    void Build3DVertex_Adaptive();

    int run;
    int event;

    std::vector<Line*> *Collection_Tracks_xz;
    std::vector<Line*> *Collection_Tracks_yz;
    std::vector<Vertex3D*> *v_vertex;
    std::vector<Vertex2D*> *v_v2d_xz;
    std::vector<Vertex2D*> *v_v2d_yz;

    int ntrdtrack_x,ntrdtrack_y;
    int nvertex_2d_x;
    int nvertex_2d_y;
    int nvertex_3d;
    double vertex_x,vertex_y,vertex_z;
    double vertex_x_err,vertex_y_err,vertex_z_err;
    double vertex_chi2,vertex_chi2_x,vertex_chi2_y;
    int vertex_ntrack,vertex_ntrack_x,vertex_ntrack_y;
    int vertex_nhit,vertex_nhit_x,vertex_nhit_y;
    int vertex_is2d;

    TCanvas* c;
    void ResetVariable();
    void MergeTracks(std::vector<Line*> *v_track);
};

#endif
