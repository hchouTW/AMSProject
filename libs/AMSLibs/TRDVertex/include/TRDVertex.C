#include "TRDVertex.h"
#include "root.h"

using namespace std;
bool TRDVertex::DEBUG=0;

int Layer_Z[20]={86,89,92,95,
                 97,100,103,106,
                 109,112,115,118,
                 121,124,127,130,
                 132,135,138,142};   // Just for fun



bool comp2DVertexByNTracks(Vertex2D *a, Vertex2D *b){
    return a->GetNTracks()>b->GetNTracks();
}

bool comp3DVertexByNTracks(Vertex3D *a, Vertex3D *b){

    if(a->GetNTracks()!=b->GetNTracks())return a->GetNTracks()>b->GetNTracks();
    else return (a->GetNhitX()+a->GetNhitY()) > (b->GetNhitX()+b->GetNhitY());
}

bool comp2DVertexByNHits(Vertex2D *a, Vertex2D *b){
    return a->GetTotalNHits()>b->GetTotalNHits();
}


bool LessThanThreeHits(Line* line){
    return line->GetNhit()<=3;
}

bool compByState(Cell* a, Cell* b){
    return a->state>b->state;
}

bool compByLength(const Cell *a, const Cell *b){
    return a->length>b->length;
}



void TRDVertex::SetRandomHits(){

    //    for(int ilayer=0;ilayer<20;ilayer++){

    //        int nhit=gRandom->Uniform(1,3);

    //        for(int i=0;i<nhit;i++){
    //            Hits.Add(TRD2DHit(gRandom->Uniform(0,20),Layer_Z[ilayer],100,ilayer));
    //        }

    //    }

    //    Hits.Add(TRD2DHit(0,Layer_Z[4],100,4));
    //    Hits.Add(TRD2DHit(1,Layer_Z[5],100,5));
    //    Hits.Add(TRD2DHit(4,Layer_Z[6],100,6));
    //    Hits.Add(TRD2DHit(5,Layer_Z[7],100,7));
    //    Hits.Add(TRD2DHit(9,Layer_Z[8],100,8));
    //    Hits.Add(TRD2DHit(12,Layer_Z[9],100,9));
    //    Hits.Add(TRD2DHit(16,Layer_Z[10],100,10));
    //    Hits.Add(TRD2DHit(17,Layer_Z[11],100,11));
    //    Hits.Add(TRD2DHit(10,Layer_Z[4],100,4));
    //    Hits.Add(TRD2DHit(14,Layer_Z[5],100,5));
    //    Hits.Add(TRD2DHit(17,Layer_Z[6],100,6));
    //    Hits.Add(TRD2DHit(19,Layer_Z[7],100,7));
    //    Hits.Add(TRD2DHit(24,Layer_Z[8],100,8));
    //    Hits.Add(TRD2DHit(26,Layer_Z[9],100,9));
    //    Hits.Add(TRD2DHit(29,Layer_Z[10],100,10));
    //    Hits.Add(TRD2DHit(32,Layer_Z[11],100,11));

    //    Hits.Add(TRD2DHit(-80,Layer_Z[4],100,4));
    //    Hits.Add(TRD2DHit(-12,Layer_Z[5],100,5));
    //    Hits.Add(TRD2DHit(4,Layer_Z[6],100,6));
    //    Hits.Add(TRD2DHit(-5,Layer_Z[7],100,7));

    //    Hits.Add(TRD2DHit(-30,Layer_Z[4],100,4));
    //    Hits.Add(TRD2DHit(-82,Layer_Z[5],100,5));
    //    Hits.Add(TRD2DHit(64,Layer_Z[6],100,6));
    //    Hits.Add(TRD2DHit(-35,Layer_Z[7],100,7));
    //    Hits.Add(TRD2DHit(-2,Layer_Z[8],100,8));
    //    Hits.Add(TRD2DHit(4,Layer_Z[9],100,9));
    //    Hits.Add(TRD2DHit(-24,Layer_Z[10],100,10));






    //    Hits.Add(TRD2DHit(3,Layer_Z[0],100,0));
    //    Hits.Add(TRD2DHit(3,Layer_Z[1],100,1));
    //    Hits.Add(TRD2DHit(3,Layer_Z[2],100,2));
    //    Hits.Add(TRD2DHit(4,Layer_Z[3],100,3));

    //    Hits.Add(TRD2DHit(7,Layer_Z[0],100,0));
    //    Hits.Add(TRD2DHit(7,Layer_Z[1],100,1));
    //    Hits.Add(TRD2DHit(7,Layer_Z[2],100,2));
    //    Hits.Add(TRD2DHit(6,Layer_Z[3],100,3));

    //    Hits.Add(TRD2DHit(5,Layer_Z[2],100,2));





    //    Hits.Sort();
}


void TRDVertex::SetHit(std::vector<SimpleHitState> *MySimpleHits_TRD){

    if(Hits.NTotalHits!=0)ClearMemory();

    for(int i=0;i<MySimpleHits_TRD->size();i++){
        //        TRD2DHit hit(MySimpleHits_TRD->at(i));
        //            hit.Print();
        //            Hits.Add(MySimpleHits_TRD->at(i));
        Hits.Add_WidthMerging(MySimpleHits_TRD->at(i));
    }
    Hits.Sort();
}

void TRDVertex::Init()
{


    Collection_Tracks_xz = new vector<Line*>;
    Collection_Tracks_yz = new vector<Line*>;
    v_vertex = new vector<Vertex3D*>;
    v_v2d_xz = new vector<Vertex2D*>;
    v_v2d_yz = new vector<Vertex2D*>;

    DEBUG=0;
    run=0;
    event=0;
    ntrdtrack_x=0;
    ntrdtrack_y=0;
    nvertex_2d_x=0;
    nvertex_2d_y=0;
    nvertex_3d=0;
    vertex_x=0;
    vertex_y=0;
    vertex_z=0;
    vertex_x_err=0;
    vertex_y_err=0;
    vertex_z_err=0;
    vertex_chi2=0;
    vertex_chi2_x=0;
    vertex_chi2_y=0;
    vertex_ntrack=0;
    vertex_ntrack_x=0;
    vertex_ntrack_y=0;
    vertex_nhit=0;
    vertex_nhit_x=0;
    vertex_nhit_y=0;
    vertex_is2d=0;




}

void TRDVertex::ResetVariable()
{
    run=0;
    event=0;
    ntrdtrack_x=0;
    ntrdtrack_y=0;
    nvertex_2d_x=0;
    nvertex_2d_y=0;
    nvertex_3d=0;
    vertex_x=0;
    vertex_y=0;
    vertex_z=0;
    vertex_x_err=0;
    vertex_y_err=0;
    vertex_z_err=0;
    vertex_chi2=0;
    vertex_chi2_x=0;
    vertex_chi2_y=0;
    vertex_ntrack=0;
    vertex_ntrack_x=0;
    vertex_ntrack_y=0;
    vertex_nhit=0;
    vertex_nhit_x=0;
    vertex_nhit_y=0;
    vertex_is2d=0;
}


int TRDVertex::Reconstruction(AMSEventR * pev){
    if(!pev) return 0;
    for(int i=0;i<pev->NTrdRawHit();i++){
        if(!pev->pTrdRawHit(i))continue;
        TRD2DHit *hit = new TRD2DHit(pev->pTrdRawHit(i));
        Hits.Add_WidthMerging(hit);
    }
    Hits.Sort();

    return  Reconstruction();
}


int TRDVertex::Reconstruction()
{
    ResetVariable();

    if(!Hits.NTotalHits)return 0;
    TrackFinding();
    ntrdtrack_x=Collection_Tracks_xz->size();
    ntrdtrack_y=Collection_Tracks_yz->size();
    VertexFinding_New();
    ClearMemory();
    return 1;
}



void TRDVertex::MakeCell_New_Improved(std::vector<Cell*>& cells,TString Projection)
{


    static int layer_y[8]={0,1,2,3,16,17,18,19};
    static int layer_x[12]={4,5,6,7,8,9,10,11,12,13,14,15};
    static int reach=5;



    int nlay;
    int layer[12];
    double x_limit=50;
    if(Projection.Contains("x") ||Projection.Contains("X")){
        nlay=12;
        for(int i=0;i<nlay;i++)layer[i]=layer_x[i];
    }else{
        nlay=8;
        for(int i=0;i<nlay;i++)layer[i]=layer_y[i];
    }

    std::vector<Cell*> v_candidates;

    for(int layer1=0;layer1<nlay-1;layer1++){
        for(int i1=0;i1<Hits.Size(layer[layer1]);i1++){
            TRD2DHit* hit1=Hits.GetHit(layer[layer1],i1);
            for(int layer2=layer1;layer2<layer1+reach+1;layer2++){
                if(layer2>nlay-1)continue;
                for(int i2=0;i2<Hits.Size(layer[layer2]);i2++){
                    TRD2DHit* hit2=Hits.GetHit(layer[layer2],i2);
                    if(hit1==hit2)continue;
                    if(hit1->z>=hit2->z)continue;
                    if(hit1->GetDistance(hit2)>100)continue;
                    Cell* cnew=new Cell(hit1,hit2);
                    v_candidates.push_back(cnew);
                    //                    hit1->AddConnectedHitUpper(hit2);
                    hit1->AddConnectedUpper(cnew);
                    //                    hit2->AddConnectedHitLower(hit1);
                    hit2->AddConnectedLower(cnew);
                }
            }
        }
    }
    if(DEBUG)cout<<"v_candidates size: " <<v_candidates.size()<<endl;;




    //======Keep Shortest connection=========
    for(int layer1=0;layer1<nlay-1;layer1++){
        for(int i1=0;i1<Hits.Size(layer[layer1]);i1++){
            TRD2DHit* hit1=Hits.GetHit(layer[layer1],i1);

            double min_length=100;
            Cell* cmin=NULL;
            for(int icell=0;icell<hit1->v_connected_upper.size();icell++){
                Cell* c=hit1->v_connected_upper.at(icell);
                if(c->length<min_length){
                    min_length=c->length;
                    cmin=c;
                }
            }

            if(cmin){
                cells.push_back(cmin);
            }

        }
    }

    for(int layer1=0;layer1<nlay-1;layer1++){
        for(int i1=0;i1<Hits.Size(layer[layer1]);i1++){
            TRD2DHit* hit1=Hits.GetHit(layer[layer1],i1);

            double min_length=100;
            Cell* cmin=NULL;
            for(int icell=0;icell<hit1->v_connected_lower.size();icell++){
                Cell* c=hit1->v_connected_lower.at(icell);
                if(c->length<min_length){
                    min_length=c->length;
                    cmin=c;
                }
            }

            if(cmin){
                if(std::find(cells.begin(),cells.end(),cmin)==cells.end()){
                    cells.push_back(cmin);
                }
            }
        }
    }

    //======Connect best matching het for every cell=========

    vector<Cell*> cells_add;

    for(double cos_cut=0.99;cos_cut>0.9;cos_cut-=0.01){
        for(int icell=0;icell<cells.size();icell++){
            Cell* c=cells.at(icell);

            TRD2DHit* p1=(TRD2DHit*)c->p1;

            if(p1){
                double min_length=100;
                Cell* cbest=NULL;
                for(int i=0;i<p1->v_connected_lower.size();i++){
                    Cell* c2=p1->v_connected_lower.at(i);
                    double costheta=c->CosTheta(c2);
                    if(costheta>cos_cut && c2->length < min_length){
                        min_length = c2->length;
                        cbest=c2;
                    }
                }
                if(cbest && std::find(cells_add.begin(),cells_add.end(),cbest)==cells_add.end()){
                    cells_add.push_back(cbest);
                }
            }

            TRD2DHit* p2=(TRD2DHit*)c->p2;

            if(p2){
                double min_length=100;
                Cell* cbest=NULL;
                for(int i=0;i<p2->v_connected_upper.size();i++){
                    Cell* c2=p2->v_connected_upper.at(i);
                    double costheta=c->CosTheta(c2);
                    if(costheta>cos_cut && c2->length < min_length){
                        min_length = c2->length;
                        cbest=c2;
                    }
                }
                if(cbest && std::find(cells_add.begin(),cells_add.end(),cbest)==cells_add.end()){
                    cells_add.push_back(cbest);
                }
            }
        }

        for(int i=0;i<cells_add.size();i++){
            if(std::find(cells.begin(),cells.end(),cells_add.at(i))==cells.end()){
                cells.push_back(cells_add.at(i));
            }
        }


    }
    //    cout<<"Cell_Add: "<<cells_add.size()<<endl;

    for(int i=0;i<v_candidates.size();i++){
        Cell* cell=(v_candidates.at(i));
        if(std::find(cells.begin(),cells.end(),cell)==cells.end())delete cell;
    }



    if(DEBUG)cout<<"Cells: " <<cells.size()<<endl;

    for(int layer1=0;layer1<nlay-1;layer1++){
        for(int i1=0;i1<Hits.Size(layer[layer1]);i1++){
            TRD2DHit* hit1=Hits.GetHit(layer[layer1],i1);
            hit1->ClearConnections();
        }
    }
    for(int i=0;i<cells.size();i++){
        Cell* cell=(cells.at(i));
        TRD2DHit* hit1=(TRD2DHit*)cell->p1;
        TRD2DHit* hit2=(TRD2DHit*)cell->p2;
        hit1->AddConnectedUpper(cell);
        hit2->AddConnectedLower(cell);
    }
    for(int i=0;i<cells.size();i++){
        Cell* c=(cells.at(i));
        c->AddConnectedCell();
    }


    return;


}


void TRDVertex::MakeCell_New(std::vector<Cell*>& cells,TString Projection)
{

    static int layer_y[8]={0,1,2,3,16,17,18,19};
    static int layer_x[12]={4,5,6,7,8,9,10,11,12,13,14,15};
    static int reach=2;



    int nlay;
    int layer[12];
    double x_limit=50;
    if(Projection.Contains("x") ||Projection.Contains("X")){
        nlay=12;
        for(int i=0;i<nlay;i++)layer[i]=layer_x[i];
    }else{
        nlay=8;
        for(int i=0;i<nlay;i++)layer[i]=layer_y[i];
    }

    vector<Cell*> v_candidates;
    for(int layer1=0;layer1<nlay-1;layer1++){
        for(int i1=0;i1<Hits.Size(layer[layer1]);i1++){
            TRD2DHit* hit1=Hits.GetHit(layer[layer1],i1);
            for(int layer2=layer1;layer2<layer1+reach+1;layer2++){
                if(layer2>nlay-1)continue;
                //                if(abs(layer[layer2]-layer[layer1])>5 )continue;
                for(int i2=0;i2<Hits.Size(layer[layer2]);i2++){
                    TRD2DHit* hit2=Hits.GetHit(layer[layer2],i2);
                    if(hit1==hit2)continue;
                    if(hit1->z>=hit2->z)continue;
                    //                    if(hit1->layer==hit2->layer && hit1->GetDistance(hit2)>5)continue;
                    v_candidates.push_back(new Cell(hit1,hit2));
                }
            }
        }
    }
    if(DEBUG)cout<<"v_candidates size: " <<v_candidates.size()<<endl;;


    //======== Only Keep best incomming and outgoing connection
    vector<Cell*> v_candidates_reduced;

    //    cout<<"Number of Connections: " <<v_candidates.size()<<endl;

    vector<Cell*>::iterator iter1=v_candidates.begin();
    vector<Cell*>::iterator iter2=v_candidates.begin();

    vector<Cell*>::iterator iter_best_incomming;
    vector<Cell*>::iterator iter_best_outgoing;


    int v_candidate_size=v_candidates.size();
    std::sort(v_candidates.begin(),v_candidates.end(),compByLength);
    std::reverse(v_candidates.begin(),v_candidates.end());



    double T=1;
    double Tc=0;
    int count=0;
    float dummy_x,dummy_y;

    while (T<100){
        count++;
        Tc=0;
        for( iter1=v_candidates.begin();iter1!=v_candidates.end();){
            Cell* c1=(*iter1);
            int iscrossed=0;
            int hassamesatrt=0;
            if(c1->length>T)break;
            for(int j=0;j<v_candidates_reduced.size();j++){
                Cell* c2=(v_candidates_reduced.at(j));
                if(!c2->IsValid())continue;
                if(c1==c2)continue;
                if(c1->p1==c2->p1  && c1->p2->layer==c2->p2->layer && c1->length>2*c2->length){
                    hassamesatrt=1;
                    break;
                }
                if(c1->Intersect(c2,dummy_x,dummy_y)){
                    iscrossed=1;
                    break;
                }
            }
            if(!iscrossed && !hassamesatrt){
                v_candidates_reduced.push_back(c1);
                v_candidates.erase(iter1);
            }else{
                iter1++;
            }
        }
        T*=2;
    }
    //    v_candidates_reduced=v_candidates;
    if(DEBUG)cout<<"Number of Reduced Cells:  " <<v_candidates_reduced.size()<<endl;


    for( iter1=v_candidates_reduced.begin();iter1!=v_candidates_reduced.end();iter1++){
        double best_weight_incomming=0;
        double best_weight_outgoing=0;
        int n_incomming=0;
        int n_outgoing=0;
        for( iter2=v_candidates_reduced.begin();iter2!=v_candidates_reduced.end();iter2++){
            if((*iter1)->p1==(*iter2)->p2){
                n_incomming++;
                double weight=(*iter1)->GetTransformFunction((*iter2),0,1);
                if(weight>best_weight_incomming){
                    best_weight_incomming=weight;
                    iter_best_incomming=iter2;
                }
            }else if((*iter1)->p2==(*iter2)->p1){
                n_outgoing++;
                double weight=(*iter1)->GetTransformFunction((*iter2),0,1);
                if(weight>best_weight_outgoing){
                    best_weight_outgoing=weight;
                    iter_best_outgoing=iter2;
                }
            }
        }
        if(!n_incomming || !n_outgoing)(*iter1)->state++;  // cell with no incomming or outgoing, increase my self
        if(best_weight_incomming)(*iter_best_incomming)->state++;
        if(best_weight_outgoing)(*iter_best_outgoing)->state++;
    }

    for(int i=0;i<v_candidates_reduced.size();i++){
        if(v_candidates_reduced.at(i)->state>1){
            v_candidates_reduced.at(i)->state=1;
            cells.push_back(v_candidates_reduced.at(i));
        }else{
            delete (v_candidates_reduced.at(i));
        }
    }


    if(DEBUG)cout<<"Final Cells: " <<cells.size()<<endl;


    for(int i=0;i<cells.size();i++){
        Cell* cell=cells.at(i);
        TRD2DHit* hit1=cell->p1;
        TRD2DHit* hit2=cell->p2;
        hit1->AddConnectedUpper(cell);
        hit2->AddConnectedLower(cell);
    }



    //=======Add Connected Cell
    for(int i=0;i<cells.size();i++){
        Cell* c=cells.at(i);
        c->AddConnectedCell();
        //        c->Print();
    }

    if(DEBUG)cout<<"NCells: " <<cells.size()<<endl;
}

void TRDVertex::Evolve_NN(std::vector<Cell*> &cells,int iside)
{

    if(DEBUG)cout<<"Evolve NN.............."<<endl;

    int cell_size=cells.size();

    gRandom->SetSeed(run+100*event);

    double Temperature=1000;
    double gamma=20;
    double beta=-0.5;
    double alpha=1;

    for(int i=0;i<cell_size;i++){
        Cell* cell_i=cells[i];
        cell_i->state=gRandom->Uniform(0.2,0.8);
    }

    double Tij,fi,fj;

    Cell* cell_i=0;
    Cell* cell_j=0;

    int sweep=0;
    double delta=100;
    double scale=0.5;

    int repeat=0;
    int Niter=0;
    while( delta>1e-5  || repeat<5){
        if(Niter>1000)break;
        Niter++;

        if(delta<1e-5){
            repeat++;
            alpha*=2;
        }else{
            Temperature*=scale;
        }
        //        beta*=scale;
        sweep++;
        //        double total_state=0;
        //        if(sweep>1) g_Energy->SetPoint(g_Energy->GetN(),sweep,Energy    );
        //        for(int i=0;i<cell_size;i++){
        //  for(cell_i=cells.begin();cell_i!=cells.end();cell_i++){
        //            Cell* cell_i=&(cells[i]);
        //            if(i<100) g_V[i]->SetPoint(g_V[i]->GetN(),sweep,cell_i->state);
        //        if(cell_i->state>0)total_state+=1;
        //        }
        //        cout<<"================================"<<endl;
        //        cout<<"Iter: "<<sweep<<endl;
        //        cout<<"Temperature: " <<Temperature<<endl;
        //        cout<<"T: " <<T<<" ,  Energy: " <<Energy<<endl;
        //        cout<<"Total State: " <<total_state<<",  Nhit: " <<Nhits<<endl;

        //Fast track finding with neural networks
        //Georg Stimpfl-Abele  1991
        //Computer Physics Communications
        delta=0;

        //        double current_Energy=0;
        for(int i=0;i<cell_size;i++){
            cell_i=cells[i];

            //            if(fabs(cell_i->p1->x-23)<1 && fabs(cell_i->p1->z-97)<1 && fabs(cell_i->p2->x-(-2.2))<0.2 && fabs(cell_i->p2->z-98.3)<0.2){
            //                cell_i->Print();
            //                cout<<"Conection : "<<endl;
            //            }



            //            DoPrint=0;
            //            if(DoPrint){cout<<"cell i : ";cell_i->Print();}
            cell_i->evolve_next_generation=0;
            fi=cell_i->state;
            for(int j=0;j<cell_size;j++){
                if(i==j)continue;
                cell_j=cells[j];
                Tij=cell_i->GetTransformFunction(cell_j,alpha,gamma);
                fj=cell_j->state;

                cell_i->evolve_next_generation+=Tij*fj;

                //                if(fabs(cell_i->p1->x-23)<1 && fabs(cell_i->p1->z-97)<1 && fabs(cell_i->p2->x-(-2.2))<0.3 && fabs(cell_i->p2->z-98.3)<0.3){
                //                    if(Tij!=0){
                //                        cell_j->Print();
                //                        cout<<"CosTheta: " <<cell_i->theta<<", "<<cell_j->theta<<",   "<<cell_i->GetAbsoluteAngle(cell_j)<<",  "<<fabs(cos(cell_i->GetAbsoluteAngle(cell_j)))<<endl;
                //                        cout<<"Tij: " <<Tij<<",   "<< cell_i->evolve_next_generation<<endl;
                //                    }

                //                }

                //                current_Energy+=-0.5*Tij*fi*fj;
            }
            cell_i->evolve_next_generation+=beta;
            //            if(DoPrint)cout<<"Change to  : "<<cell_i->evolve_next_generation<<endl;
            double newstate= 0.5*(1+TMath::TanH(cell_i->evolve_next_generation/Temperature));
            delta+=fabs(newstate-cell_i->state);


            cell_i->state =newstate ;


            //            if(fabs(cell_i->p1->x-23)<1 && fabs(cell_i->p1->z-97)<1 && fabs(cell_i->p2->x-(-2.2))<0.2 && fabs(cell_i->p2->z-98.3)<0.2){
            //                cout<<"New state: " <<endl;
            //                cell_i->Print();
            //            }


        }
        delta/=cell_size;

        //        current_Energy+=beta;
        //        dE=Energy-current_Energy;
        //        Energy=current_Energy;

        //        if(dE<0){
        //            cout<<"Energy Reversed..............Abort"<<endl;
        //            break;
        //        }



        //        for(int i=0;i<cell_size;i++){
        //            Cell* cell_i=&(cells.at(i));
        //            cell_i->state =  0.5*(1+TMath::TanH(cell_i->evolve_next_generation/Temperature));


        //            //            if(cell_i->state<0)cell_i->state=0;


        //            cell_i->evolve_next_generation=0;
        //        }




        //        TCanvas *ctest = new TCanvas(Form("c_%i",iter),Form("c_%i",iter),800,800);
        //        ctest->DrawFrame(-100,-150,100,160);
        //        for(int i=0;i<cell_size;i++){
        //            Cell* cell_i=&(cells.at(i));
        //            if(cell_i->state>0.1)cell_i->Draw();
        //        }

    }


    std::vector<Cell*>::iterator iter;
    for(iter=cells.begin();iter!=cells.end();iter++){
        if((*iter)->state > 0.5){
            (*iter)->state=1;
            //            continue;
        }else{
            (*iter)->state=-1;
            (*iter)->InValidate();
            //            cells.erase(iter);
        }
    }

    return;

}



void TRDVertex::DrawCells(std::vector<Cell*> &cells)
{
    for(int i=0;i<cells.size();i++){
        Cell* cell1=cells.at(i);
        cell1->Draw();
    }
}


void TRDVertex::Build3DVertex_Adaptive()
{

    if(DEBUG){
        cout<<"===================================="<<endl;
        cout<<"=== Build3D Vertex Adaptive :======="<<endl;
        cout<<"===================================="<<endl;
    }

    vector<Line*>::iterator it_vtr;

    std::vector<Line*> _Collection_Tracks_xz = *Collection_Tracks_xz;
    std::vector<Line*> _Collection_Tracks_yz = *Collection_Tracks_yz;

    while(_Collection_Tracks_xz.size()>1 ||  _Collection_Tracks_yz.size()>1){
        if(DEBUG)cout<<"XY N Track candidate: " <<_Collection_Tracks_xz.size()<<",  "<<_Collection_Tracks_yz.size()<<endl;
        Vertex3D* candidate = new Vertex3D;
        if(_Collection_Tracks_xz.size()>=2){
            for(int i=0;i<_Collection_Tracks_xz.size();i++){
                candidate->AddTrack_XZ(_Collection_Tracks_xz.at(i));
            }
        }
        if(_Collection_Tracks_yz.size()>=2){
            for(int i=0;i<_Collection_Tracks_yz.size();i++){
                candidate->AddTrack_YZ(_Collection_Tracks_yz.at(i));
            }
        }
        if(DEBUG)cout<<"N Track : "<<candidate->v_x->GetNTracks()<<",   "<<candidate->v_y->GetNTracks()<<endl;
        int result=candidate->CombineFit_3D_Adaptive();
        if(result>0){
            v_vertex->push_back(candidate);
            if(DEBUG){
                cout<<"Candidate: "<<v_vertex->size()<<endl;
                candidate->Print();
            }

            for(it_vtr=_Collection_Tracks_xz.begin();it_vtr!=_Collection_Tracks_xz.end();){
                if(candidate->v_x->Contains(*it_vtr)){
                    _Collection_Tracks_xz.erase(it_vtr);
                }
                else it_vtr++;
            }

            for(it_vtr=_Collection_Tracks_yz.begin();it_vtr!=_Collection_Tracks_yz.end();){
                if(candidate->v_y->Contains(*it_vtr)){
                    _Collection_Tracks_yz.erase(it_vtr);
                }else it_vtr++;
            }

        }else{
            delete candidate;
            break;
        }
    }

    if(DEBUG){
        cout<<"=========================================="<<endl;
        cout<<"3D Vertex Candidate: " <<v_vertex->size()<<endl;
        for(int i=0;i<v_vertex->size();i++){
            v_vertex->at(i)->Print();
        }
        cout<<"=========================================="<<endl;
    }
    return;


}


void TRDVertex::Build2DVertex_Adaptive(int iside)
{

    if(DEBUG){
        cout<<"===================================="<<endl;
        cout<<"===Build2DVertex_Adaptive : "<<iside<<"========"<<endl;
        cout<<"===================================="<<endl;

    }


    vector<Line*> v_track;
    vector<Vertex2D*> *v_v2d;

    if(iside==0){
        v_v2d=v_v2d_xz;
        v_track=*Collection_Tracks_xz;
    }else{
        v_v2d=v_v2d_yz;
        v_track=*Collection_Tracks_yz;
    }
    vector<Line*>::iterator it_vtr;

    while(v_track.size()>1){
        cout<<"N Track candidate: " <<v_track.size()<<endl;
        Vertex2D* candidate = new Vertex2D;
        for(int i=0;i<v_track.size();i++){
            candidate->Add(v_track.at(i));
        }
        candidate->CombineFit_New_Adaptive();

        cout<<"Candidate: "<<v_v2d->size()<<endl;
        candidate->Print();

        v_v2d->push_back(candidate);

        for(it_vtr=v_track.begin();it_vtr!=v_track.end();){
            if(candidate->Contains(*it_vtr))v_track.erase(it_vtr);
            else it_vtr++;
        }
    }

    if(DEBUG)cout<<"2D Vertex Candidate: " <<v_v2d->size()<<endl;

    return;
}

void TRDVertex::Build2DVertex_ZVTOP(int iside)
{

    vector<Line*> *v_track;
    vector<Vertex2D*> *v_v2d;

    if(iside==0){
        v_v2d=v_v2d_xz;
        v_track=Collection_Tracks_xz;
    }else{
        v_v2d=v_v2d_yz;
        v_track=Collection_Tracks_yz;
    }

    VertexFinder2D_ZVTOP *finder=new VertexFinder2D_ZVTOP(iside);
    finder->SetLineColection(v_track);
    finder->SetVertexColection(v_v2d);
    finder->FindVertex();

    if(DEBUG)cout<<"2D Vertex Candidate: " <<v_v2d->size()<<endl;


    delete finder;
    return;
}

void TRDVertex::Evolve(std::vector<Cell*> &cells,int iside)
{
    int maxNiter=500;


    for(int iter=0;iter<maxNiter;iter++){
        int flag_changed=0;
        for(int i=0;i<cells.size();i++){
            Cell* cell1=(cells.at(i));
            cell1->evolve_next_generation=0;
            if(cell1->state<0)continue;
            if(!cell1->v_connected_upper.size())continue;
            for(int j=0;j<cell1->v_connected_upper.size();j++){
                Cell* cell2=cell1->v_connected_upper.at(j);
                if(cell2->state<0)continue;
                if(cell1->IsSameState(cell2)){
                    cell1->evolve_next_generation=1;
                    flag_changed=1;
                    break;
                }
            }
        }

        if(flag_changed){
            for(int i=0;i<cells.size();i++){
                Cell* cell1=(cells.at(i));
                cell1->Evolve();
            }
        }else {
            break;
        }

    }

    return;
}

void TRDVertex::MergeTracks(vector<Line*> *v_track)
{
    if(DEBUG)cout<<"Total Tracks Segment: "<<v_track->size()<<endl;
    int maxiter=10;
    int iter=0;
    int oldsize=v_track->size();

    vector<Line*>::iterator it_vtr;
    vector<Line*>::iterator it2_vtr;
    while(iter<maxiter && v_track->size()){
        for(it_vtr=v_track->begin();it_vtr!=v_track->end();it_vtr++){
            for(it2_vtr=it_vtr+1;it2_vtr!=v_track->end();){
                if((*it_vtr)->Merge(*it2_vtr)){
                    delete (*it2_vtr);
                    v_track->erase(it2_vtr);
                }
                else it2_vtr++;
            }
        }
        if(v_track->size()==oldsize){
            break;
        }
        oldsize=v_track->size();
        iter++;
    }

    // v_track->erase( std::remove_if(v_track->begin(), v_track->end(), LessThanThreeHits), v_track->end() );
    //    vector<Line*>::iterator it;
    //    for(it=v_track->begin();it!=v_track->end();){
    //        if(LessThanThreeHits(*it)){
    //            delete *it;
    //            it=v_track->erase(it);
    //        }else{
    //            it++;
    //        }
    //    }

    if(DEBUG){
        cout<<"Total Tracks: "<<v_track->size()<<endl;
        for(int i=0;i<v_track->size();i++){
            v_track->at(i)->Draw();
            v_track->at(i)->Print();
        }
    }
}

void TRDVertex::SelectTracks(std::vector<Cell*>& cells,int iside){


    // Need to sort the pointer,  not the object, otherwise the point ill be messed up !!!
    vector<Cell*> v_p_cells_sorted=cells;
    //    for(int i=0;i<cells.size();i++)v_p_cells_sorted.push_back((cells.at(i)));   //Make a vector of pointer pointing to the original cells,
    std::sort(v_p_cells_sorted.begin(),v_p_cells_sorted.end(),compByState);

    //=======Init for track and vertex

    vector<Line*> *v_track;
    vector<Vertex2D*> *v_v2d;

    if(iside==0){
        v_v2d=v_v2d_xz;
        v_track=Collection_Tracks_xz;
    }else{
        v_v2d=v_v2d_yz;
        v_track=Collection_Tracks_yz;
    }

    if(DEBUG)cout<<"Cells Size: " <<v_p_cells_sorted.size()<<endl;
    //    for(int i=0;i<v_p_cells_sorted.size();i++){
    //        Cell* cell=v_p_cells_sorted.at(i);
    //        cell->Print();
    //    }

    for(int i=0;i<v_p_cells_sorted.size();i++){
        //    for(int i=0;i<2;i++){
        Cell* cell=v_p_cells_sorted.at(i);
        if(cell->assigned)continue;
        if(cell->state<0)continue;
        //        cell->Print();
        Segment newsegment(cell);
        //        cout<<"Setment nhits: " <<newsegment.GetNHits()<<", Length : "<<newsegment.GetLength()<<",  nhits: " <<newsegment.GetNHits()<<endl;
        if(newsegment.GetNHits()>=3/* && newsegment.GetLength()/newsegment.GetNHits()<10*/){
            double chi2, ndof;
            newsegment.FitStraightLine(chi2,ndof);
            if(DEBUG)cout<<"Chi2: " <<chi2<<"/"<<ndof<<endl;
            //            if(chi2/ndof<100){
            //                newsegment.Print();
            //                newsegment.Draw();
            Line* par  = new Line(iside);
            for(int i=0;i<newsegment.GetNHits();i++ ) {
                par->AddHit(newsegment.hits.at(i));
            }
            int result=par->LineFit_TMinuit_Robust();
            v_track->push_back(par);
            for(int icell=0;icell<newsegment.cells.size();icell++){
                (newsegment.cells.at(icell))->InValidateOtherCells();  // This needs improvement  maybe after robust fit ?
            }

            //            }
        }
    }

    //    return;
    //    if(v_track->size()<2)return;

    //    DEBUG=1;

    //=====Combine different track segment into one single track==========

    MergeTracks(v_track);




}




void TRDVertex::SelectTracks_New(std::vector<Cell*>& cells,int iside){


    vector<Cell*> v_p_cells_sorted=cells;
    std::sort(v_p_cells_sorted.begin(),v_p_cells_sorted.end(),compByState);

    //=======Init for track and vertex

    vector<Line*> *v_track;
    vector<Vertex2D*> *v_v2d;

    if(iside==0){
        v_v2d=v_v2d_xz;
        v_track=Collection_Tracks_xz;
    }else{
        v_v2d=v_v2d_yz;
        v_track=Collection_Tracks_yz;
    }

    if(DEBUG)cout<<"Cells Size: " <<v_p_cells_sorted.size()<<endl;


    for(double length_threshold=12;length_threshold>=3;length_threshold--){
        if(DEBUG)cout<<"Length : " <<length_threshold<<endl;

        double costheta_cut=0.1;

        while(1-costheta_cut>0.8){
            if(DEBUG)cout<<"====Cos Theta: " <<1-costheta_cut<<endl;
            for(int i=0;i<v_p_cells_sorted.size();i++){
                Cell* cell=v_p_cells_sorted.at(i);
                if(cell->assigned)continue;
                if(cell->state<0)continue;
                Segment newsegment(cell,1-costheta_cut);
                if(newsegment.GetNHits()>=length_threshold){
                    double chi2, ndof;
                    newsegment.FitStraightLine(chi2,ndof);
                    if(DEBUG)cout<<"Chi2: " <<chi2<<"/"<<ndof<<endl;
                    Line* par  = new Line(iside);
                    for(int i=0;i<newsegment.GetNHits();i++ ) {
                        par->AddHit(newsegment.hits.at(i));
                    }
                    int result=par->LineFit_TMinuit_Robust();

                    if(par->GetTotalChi2()/(par->GetNhit()-2) <5){
                        v_track->push_back(par);
                        newsegment.Finish();
                        for(int icell=0;icell<newsegment.cells.size();icell++){
                            (newsegment.cells.at(icell))->InValidateOtherCells();  // This needs improvement  maybe after robust fit ?
                        }
                    }else{
                        delete par;
                    }
                }
            }
            costheta_cut*=2;
        }
    }

    //    cout<<"NTracks Before Merging: " <<v_track->size()<<endl;

    MergeTracks(v_track);

}


//int TRDVertex::Build3DVertex(){

//    if(DEBUG)cout<<"Build  3D Vertex:  "<<v_v2d_xz->size()<<",   "<<v_v2d_yz->size()<<endl;

//    if(!v_v2d_xz->size() && !v_v2d_yz->size())return 0;


//    //===if ntrk==2 and no match in other side , remove

//    //    if(0){
//    //        vector<Vertex2D>::iterator iter_yz_matched;


//    //        double dist_cut=5;

//    //        for(iter_xz=v_v2d_xz->begin();iter_xz!=v_v2d_xz->end();){
//    //            if(iter_xz->GetNTracks()>2){
//    //                iter_xz++;
//    //            }else{
//    //                double mindist=999;
//    //                for(iter_yz=v_v2d_yz->begin();iter_yz!=v_v2d_yz->end();iter_yz++){
//    //                    {
//    //                        double dist=fabs(iter_xz->Z()-iter_yz->Z());
//    //                        if(mindist>dist)mindist=dist;
//    //                    }
//    //                }
//    //                if(mindist>dist_cut)v_v2d_xz->erase(iter_xz);
//    //                else iter_xz++;
//    //            }
//    //        }

//    //        for(iter_yz=v_v2d_yz->begin();iter_yz!=v_v2d_yz->end();){
//    //            if(iter_yz->GetNTracks()>2){
//    //                iter_yz++;
//    //            }else{        double mindist=999;
//    //                for(iter_xz=v_v2d_xz->begin();iter_xz!=v_v2d_xz->end();iter_xz++){
//    //                    double dist=fabs(iter_xz->Z()-iter_yz->Z());
//    //                    if(mindist>dist)mindist=dist;
//    //                }
//    //                if(mindist>dist_cut)v_v2d_yz->erase(iter_yz);
//    //                else iter_yz++;
//    //            }
//    //        }
//    //        if(DEBUG)
//    //        {
//    //            cout<<"After 2-Track Vertex Clean Up:  "<<v_v2d_xz->size()<<",   "<<v_v2d_yz->size()<<endl;
//    //            cout<<"XZ:  "<<endl;
//    //            for(iter_xz=v_v2d_xz->begin();iter_xz!=v_v2d_xz->end();iter_xz++){
//    //                iter_xz->Print();
//    //            }
//    //            cout<<"YZ:  "<<endl;
//    //            for(iter_yz=v_v2d_yz->begin();iter_yz!=v_v2d_yz->end();iter_yz++){
//    //                iter_yz->Print();
//    //            }
//    //        }





//    //    }

//    vector<Vertex2D*>::iterator iter_xz;
//    vector<Vertex2D*>::iterator iter_yz;



//    sort(v_v2d_xz->begin(),v_v2d_xz->end(),comp2DVertexByNTracks);
//    sort(v_v2d_yz->begin(),v_v2d_yz->end(),comp2DVertexByNTracks);


//    //=========Build 3D


//    for(int i=0;i<v_v2d_xz->size();i++){
//        for(int j=0;j<v_v2d_yz->size();j++){
//            if(v_v2d_xz->at(i)->dir() * v_v2d_yz->at(j)->dir()<0) continue;
//            Vertex3D *v=new Vertex3D;
//            v->SetX(v_v2d_xz->at(i));
//            v->SetY(v_v2d_yz->at(j));
//            //            v->CombineFit();
//            //            v->RemoveOutliers();
//            v->CombineFit_New();
//            if(DEBUG){            //======Candidate
//                cout<<"=========Candidate=============="<<endl;
//                cout<<"Vertex: " <<v->x<<", "<<v->y<<", "<<v->z<<endl;
//                cout<<"Vertex err: " <<v->x_err<<", "<<v->y_err<<", "<<v->z_err<<endl;
//                cout<<"ntrack : "<<v->GetNTracksX()<<",  "<<v->GetNTracksY()<<endl;
//                cout<<"Fitted ,  chi2 : " <<v->GetChi2X()<<",   "<<v->GetChi2Y()<<endl;
//                v->Print();
//            }
//            v_vertex->push_back(v);

//            //            if(v->GetChi2X()/v->GetNhitX() <5 && v->GetChi2Y()/v->GetNhitY()<5) {
//            //                if(DEBUG){
//            //                    cout<<"3D Candidate: "<<endl;
//            //                    v->Print();
//            //                }
//            //                v_vertex->push_back(v);
//            //            }else delete v;
//        }
//    }

//    sort(v_vertex->begin(),v_vertex->end(),comp3DVertexByNTracks);



//    vector<Vertex3D*>::iterator iter_3d;
//    vector<Vertex3D*>::iterator iter2_3d;

//    for(iter_3d=v_vertex->begin();iter_3d!=v_vertex->end();iter_3d++){
//        for(iter2_3d=iter_3d+1;iter2_3d!=v_vertex->end();){
//            if((*iter_3d)->v_x == (*iter2_3d)->v_x || (*iter_3d)->v_y == (*iter2_3d)->v_y){
//                delete (*iter2_3d);
//                v_vertex->erase(iter2_3d);
//            }
//            else iter2_3d++;
//        }
//    }

//    for(iter_3d=v_vertex->begin();iter_3d!=v_vertex->end();iter_3d++){
//        for(iter2_3d=iter_3d+1;iter2_3d!=v_vertex->end();){
//            if((*iter_3d)->HasSharedTracks((*iter2_3d)) && (*iter2_3d)->GetNTracks()<=4){
//                delete (*iter2_3d);
//                v_vertex->erase(iter2_3d);
//            }            else iter2_3d++;
//        }
//    }
//    //    for(iter_xz=v_v2d_xz->begin();iter_xz!=v_v2d_xz->end();){
//    //        bool matched=0;
//    //        for(iter_yz=v_v2d_yz->begin();iter_yz!=v_v2d_yz->end();){
//    //            {
//    //                //                double dist=fabs(iter_xz->Z()-iter_yz->Z());
//    //                //                if(dist < dist_cut){

//    //                //==Remove obviously wrong direction
//    //                if((*iter_xz).dir() * (*iter_yz).dir()<0){
//    //                    iter_yz++;
//    //                    continue;
//    //                }

//    //                Vertex3D v;
//    //                v.SetX(&*(iter_xz));
//    //                v.SetY(&*(iter_yz));
//    //                v.CombineFit();
//    //                if(v.GetChi2X()/v.GetNhitX() <5 && v.GetChi2Y()/v.GetNhitY()<5){
//    //                    matched=1;
//    //                    v.Print();
//    //                    v_vertex->push_back(v);
//    //                    v_v2d_yz->erase(iter_yz);
//    //                    break;
//    //                }else{
//    //                    iter_yz++;
//    //                }
//    //            }
//    //        }
//    //        if(matched)v_v2d_xz->erase(iter_xz);
//    //        else iter_xz++;
//    //        //        DrawVertex(&(v_vertex->back()),v_vertex->size()-1);
//    //    }

//    // Remove matched 2D vertex
//    for(int i=0;i<v_vertex->size();i++){
//        for(iter_xz=v_v2d_xz->begin();iter_xz!=v_v2d_xz->end();iter_xz++){
//            if((*iter_xz) == v_vertex->at(i)->v_x){
//                //	    delete *iter_xz;
//                //	    v_v2d_xz->erase(iter_xz);
//                (*iter_xz)->has3D=1;
//            } //else iter_xz++;
//        }
//        for(iter_yz=v_v2d_yz->begin();iter_yz!=v_v2d_yz->end();iter_yz++){
//            if((*iter_yz) == v_vertex->at(i)->v_x){
//                //	    delete *iter_yz;
//                //v_v2d_yz->erase(iter_yz);
//                (*iter_yz)->has3D=1;
//            } //else iter_yz++;
//        }
//    }

//    //====Remove standalone 2_track 2d vertex


//    for(iter_xz=v_v2d_xz->begin();iter_xz!=v_v2d_xz->end();iter_xz++){
//        if((*iter_xz)->GetNTracks()<=2){
//            (*iter_xz)->isremoved=1;
//            //	delete *iter_xz;
//            //	v_v2d_xz->erase(iter_xz);
//        }//else iter_xz++;
//    }

//    for(iter_yz=v_v2d_yz->begin();iter_yz!=v_v2d_yz->end();iter_yz++){
//        if((*iter_yz)->GetNTracks()<=2){
//            (*iter_yz)->isremoved=1;
//            //	delete *iter_yz;
//            //	v_v2d_yz->erase(iter_yz);
//        }//else iter_yz++;
//    }






//    if(DEBUG){
//        cout<<"==========Final Candidate==============="<<endl;

//        cout<<"==========3D : ==============="<<endl;
//        for(int i=0;i<v_vertex->size();i++){
//            v_vertex->at(i)->Print();
//        }

//    }




//    //==========Temporary solution,  this needs improvement in the future

//    //    if(!v_vertex->size())return;
//    //    std::sort(v_vertex->begin(),v_vertex->end(),comp3DVertexByNTracks);
//    //    int primary_vertex_ntracks=v_vertex[0].GetNTracks();

//    //    if(primary_vertex_ntracks>4){
//    //        vector<Vertex3D>::iterator it_v3d;
//    //        for(it_v3d=v_vertex->begin();it_v3d!=v_vertex->end();){
//    //            if(it_v3d->GetNTracks()<=4) v_vertex->erase(it_v3d);
//    //            else it_v3d++;
//    //        }
//    //    }


//}





void TRDVertex::DrawVertex(int iside, Vertex2D *v2d, int index=0){
    TGraph* g_vertex_x=new TGraph();
    g_vertex_x->SetPoint(0,v2d->x,v2d->z);
    g_vertex_x->SetMarkerColor(kBlue);
    g_vertex_x->SetMarkerSize(3);
    g_vertex_x->SetMarkerStyle(20);
    g_vertex_x->Draw("p");
}

void TRDVertex::DrawVertex(int iside, Vertex3D *v3d, int index=0){


    if(iside==0){
        TGraph* g_vertex_x=new TGraph();
        g_vertex_x->SetPoint(0,v3d->x,v3d->z);
        g_vertex_x->SetMarkerSize(3);
        g_vertex_x->SetMarkerStyle(20);
        g_vertex_x->Draw("p");
        v3d->v_x->Draw();
    }

    if(iside==1){
        TGraph* g_vertex_y=new TGraph();
        g_vertex_y->SetPoint(0,v3d->y,v3d->z);
        g_vertex_y->SetMarkerSize(3);
        g_vertex_y->SetMarkerStyle(20);
        g_vertex_y->Draw("p");
        v3d->v_y->Draw();
    }
}


void TRDVertex::MyLabel(const char* drawstring,Double_t x,Double_t y,float size,Color_t color)
{
    TLatex *l=new TLatex(); //l.SetTextAlign(12); l.SetTextSize(tsize);
    l->SetNDC();
    l->SetTextFont(72);
    l->SetTextColor(color);
    l->SetLineWidth(2);
    l->SetTextSize(size);
    l->DrawLatex(x,y,drawstring);

}

void TRDVertex::TrackFinding(){


    MakeCell_New_Improved(cells_x,"XZ");

    TString name=Form("Event_XZ_%i_%i",run,event);

    if(DEBUG){
        c=new TCanvas(name,name,2400,1200);
        c->Divide(3,2);
        c->cd(1);
        gPad->DrawFrame(-150,-150,90,190);
        Hits.Draw_XZ();
        if(DEBUG)DrawCells(cells_x);
    }

    Evolve_NN(cells_x,0);
    if(DEBUG){
        c->cd(2);
        gPad->DrawFrame(-150,-150,90,190);
        Hits.Draw_XZ();
    }
    if(DEBUG)DrawCells(cells_x);

    Evolve(cells_x,0);
    if(DEBUG){
        c->cd(3);
        gPad->DrawFrame(-150,-150,90,190);
        Hits.Draw_XZ();
    }
    if(DEBUG)DrawCells(cells_x);

    SelectTracks_New(cells_x,0);





    // cout<<"NTracks   XZ: " <<Collection_Tracks_xz->size()<<endl;
    MakeCell_New_Improved(cells_y,"YZ");
    if(DEBUG){
        c->cd(4);
        gPad->DrawFrame(-150,-150,90,190);
        Hits.Draw_YZ();
        if(DEBUG) DrawCells(cells_y);
    }
    Evolve_NN(cells_y,1);
    if(DEBUG){
        c->cd(5);
        gPad->DrawFrame(-150,-150,90,190);
        Hits.Draw_YZ();
        if(DEBUG) DrawCells(cells_y);
    }
    Evolve(cells_y,1);
    if(DEBUG){
        c->cd(6);
        gPad->DrawFrame(-150,-150,90,190);
        Hits.Draw_YZ();
        if(DEBUG) DrawCells(cells_y);
    }
    SelectTracks_New(cells_y,1);
    // cout<<"NTracks   YZ: " <<Collection_Tracks_yz->size()<<endl;

    if(DEBUG){
        MyLabel(name,0.2,0.9,0.05,kBlue);
        c->Print("TRDTracks.pdf","pdf");
    }
    return;
}

void TRDVertex::ClearMemory()
{


    for(int i=0;i<Collection_Tracks_xz->size();i++)delete Collection_Tracks_xz->at(i);
    Collection_Tracks_xz->clear();

    for(int i=0;i<Collection_Tracks_yz->size();i++)delete Collection_Tracks_yz->at(i);
    Collection_Tracks_yz->clear();

    for(int i=0;i<v_vertex->size();i++)delete v_vertex->at(i);
    v_vertex->clear();

    for(int i=0;i<v_v2d_xz->size();i++)delete v_v2d_xz->at(i);
    v_v2d_xz->clear();

    for(int i=0;i<v_v2d_yz->size();i++)delete v_v2d_yz->at(i);
    v_v2d_yz->clear();




    //cout<<"DEBUG: " <<t++<<endl;
    for(int i=0;i<cells_x.size();i++){
        Cell* c=cells_x.at(i);
        delete c;
    }
    cells_x.clear();

    //cout<<"DEBUG: " <<t++<<endl;
    for(int i=0;i<cells_y.size();i++){
        Cell* c=cells_y.at(i);
        delete c;
    }
    cells_y.clear();


    Hits.ClearMemory();
}

void TRDVertex::VertexFinding_New(){
    Build3DVertex_Adaptive();
    //    if(DEBUG){
    //        c->cd(3);
    //        gPad->DrawFrame(-100,-150,100,160);
    //        for(int i=0;i<v_vertex->size();i++)DrawVertex(0,v_vertex->at(i),i);
    //        c->cd(4);
    //        gPad->DrawFrame(-100,-150,100,160);
    //        for(int i=0;i<v_vertex->size();i++)DrawVertex(1,v_vertex->at(i),i);
    //    }


    nvertex_2d_x=0;
    nvertex_2d_y=0;
    nvertex_3d=0;

    for(int i=0;i<v_vertex->size();i++){
        if(v_vertex->at(i)->Is2D){
            if(v_vertex->at(i)->HasXZ)nvertex_2d_x++;
            if(v_vertex->at(i)->HasYZ)nvertex_2d_y++;
        }else nvertex_3d++;

    }

    vertex_is2d=0;
    vertex_x=0;
    vertex_y=0;
    vertex_z=0;
    vertex_x_err=0;
    vertex_y_err=0;
    vertex_z_err=0;
    vertex_chi2=0;
    vertex_nhit=0;
    vertex_ntrack=0;
    vertex_chi2_x=0;
    vertex_nhit_x=0;
    vertex_ntrack_x=0;
    vertex_chi2_y=0;
    vertex_nhit_y=0;
    vertex_ntrack_y=0;

    if(nvertex_3d>0){
        Vertex3D* v=v_vertex->at(0);
        vertex_x=v->x;
        vertex_y=v->y;
        vertex_z=v->z;
        vertex_x_err=v->x_err;
        vertex_y_err=v->y_err;
        vertex_z_err=v->z_err;
        vertex_chi2=v->GetChi2X()+v->GetChi2Y();
        vertex_nhit=v->GetNhitX()+v->GetNhitY();
        vertex_ntrack=v->GetNTracks();
        vertex_chi2_x=v->GetChi2X();
        vertex_nhit_x=v->GetNhitX();
        vertex_ntrack_x=v->GetNTracksX();
        vertex_chi2_y=v->GetChi2Y();
        vertex_nhit_y=v->GetNhitY();
        vertex_ntrack_y=v->GetNTracksY();
        vertex_is2d=v->Is2D;
    }



    return;
}

void TRDVertex::VertexFinding(){


    if(Collection_Tracks_yz->size()>=3 && Collection_Tracks_xz->size()>=3){;
        Build2DVertex_ZVTOP(0);
        Build2DVertex_ZVTOP(1);
        //        Build3DVertex();


        if(DEBUG){
            c->cd(3);
            gPad->DrawFrame(-100,-150,100,160);
            for(int i=0;i<v_vertex->size();i++)DrawVertex(0,v_vertex->at(i),i);
            for(int i=0;i<v_v2d_xz->size();i++){
                if(!v_v2d_xz->at(i)->has3D && !v_v2d_xz->at(i)->isremoved)  DrawVertex(0,v_v2d_xz->at(i),i);
            }
            c->cd(4);
            gPad->DrawFrame(-100,-150,100,160);
            for(int i=0;i<v_vertex->size();i++)DrawVertex(1,v_vertex->at(i),i);
            for(int i=0;i<v_v2d_yz->size();i++){
                if(!v_v2d_yz->at(i)->has3D && !v_v2d_yz->at(i)->isremoved)            DrawVertex(1,v_v2d_yz->at(i),i);
            }
        }
        if(DEBUG){
            cout<<"3D Vertex: " <<v_vertex->size()<<endl;
            cout<<"2D Vertex: " <<v_v2d_xz->size()<<", "<<v_v2d_yz->size()<<endl;
        }

    }


    return;
}

void TRDVertex::CleanUp(std::vector<Cell*> &cells)
{

    //    cout<<"Before Clean Up: "<<cells.size()<<endl;
    vector<Cell*>::iterator it=cells.begin();
    for(;it!=cells.end();++it){
        Cell * mycell =*it;
        if(!mycell->v_connected_samestart.size())continue;
        bool is_remove=0;
        for(int i=0;i<mycell->v_connected_samestart.size();i++){
            Cell* startcell=mycell->v_connected_samestart.at(i);
            if(startcell->state<0)continue;
            if(!startcell->v_connected_upper.size())continue;
            for(int j=0;j<startcell->v_connected_upper.size();j++){
                Cell* endcell=startcell->v_connected_upper.at(j);
                if(endcell->p2 == mycell->p2){
                    is_remove=1;
                    break;
                }
            }
            if(is_remove)break;
        }
        if(is_remove){
            cout<<"Remove: ";
            //            (*it).Print();
            (*it)->state=-99;
        }
    }

    //    cout<<"After Clen Up: "<<cells.size()<<endl;
    return;
}
