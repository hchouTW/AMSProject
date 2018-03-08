#include "Cell.h"


using namespace std;



Cell::Cell(){
    evolve_next_generation=0;
    valid=1;
    p1=0;
    p2=0;
    theta=0;
    state=0;
    //    v_connected_lower.reserve(5);
    //    v_connected_upper.reserve(5);
    //    v_connected_sameend.reserve(5);
    //    v_connected_samestart.reserve(5);

    assigned=0;
    theta_range[0]=0;
    theta_range[1]=0;

    HasMatchedCells=0;
    length=0;

    connected_upper_minangle=999;
    connected_lower_minangle=999;

}

Cell::Cell(TRD2DHit *_p1, TRD2DHit *_p2):p1(_p1),p2(_p2){
    HasMatchedCells=0;
    valid=1;
    state=1;
    TRD2DHit* temp;
    if(_p2->z < _p1->z){  // Fomr here p1 should be always "below" p2
        temp=_p1;
        _p1=_p2;
        _p2=temp;
    }
    theta=TMath::ATan2((p2->y-p1->y),(p2->z-p1->z));
    if(p1->type !=p2->type){
        std::cout<<"~~~~~~~~Cell Creator~~~~~~Error ,  Hit in different projection "<<std::endl;
        cout<<"p1: " <<p1->layer<<", "<<p2->layer<<endl;
    }
    theta_range[0]=TMath::ATan2(((p2->y-0.5)-(p1->y+0.5)),(p2->z-p1->z));
    theta_range[1]=TMath::ATan2(((p2->y+0.5)-(p1->y-0.5)),(p2->z-p1->z));


    length=sqrt(pow(p1->y-p2->y,2)+pow(p1->z-p2->z,2));

    assigned=0;
    connected_upper_minangle=999;
    connected_lower_minangle=999;

    //    v_connected_lower.reserve(5);
    //    v_connected_upper.reserve(5);
    //    v_connected_sameend.reserve(5);
    //    v_connected_samestart.reserve(5);

    //v_hits.push_back(p1);
  //  v_hits.push_back(p2);

    v_connected_lower.clear();
    v_connected_upper.clear();
    v_connected_sameend.clear();
    v_connected_samestart.clear();


}

Cell::~Cell(){
//    v_connected_upper.clear();
//    v_connected_lower.clear();
//    v_connected_samestart.clear();
//    v_connected_sameend.clear();
}

bool Cell::IsConnected(Cell *cell2){
    if(p1==cell2->p2 || p2==cell2->p1)return 1;
    else return 0;
}

bool Cell::IsNeighbour(Cell *cell2){
    if(!IsSameState(cell2))return 0;
    return IsMatched(cell2);
}

bool Cell::IsSameState(Cell *cell2){
    return state==cell2->state;
}

bool Cell::IsMatched(Cell *cell2){

    return 1;


    //    // what is this ?
    //    static double threshoud=(5./180.)*TMath::Pi();
    //    double angle=GetAngle(cell2);

    //    double range[2];
    //    GetAngleRange(cell2,range);


    //    if(range[0]*range[1] < 0)  return 1;  // different sign means 0 is possible.  this happen when two hits are close to each other.


    //    if(angle+angle_err<threshoud && angle-angle_err<threshoud)return 1;

    //    range[0]=fabs(range[0]);
    //    range[1]=fabs(range[0]);

    //if(fabs(range[0])>TMath::Pi())range[0]=TMath::TwoPi()-fabs(range[0]);
    //if(fabs(range[1])>TMath::Pi())range[1]=TMath::TwoPi()-fabs(range[1]);

    //cout<<"angle: " <<angle<<endl<<",  Range: " <<range[0]<<", "<<range[1]<<endl;


    return 0;
}


//double Cell::GetAbsoluteAngle(Cell *cell2){
//    double diff=fabs(theta-cell2->theta);
//    if(diff>TMath::Pi())diff=TMath::TwoPi()-diff;
//    return diff;
//}

//double Cell::GetAngle(Cell *cell2){
//    double diff=theta-cell2->theta;
//    //    if(diff>TMath::Pi())diff=TMath::TwoPi()-diff;
//    return diff;
//}

//void Cell::GetAngleRange(Cell *cell2, double* range){
//    range[0]=theta_range[0]-cell2->theta_range[1];
//    range[1]=theta_range[1]-cell2->theta_range[0];
//}

bool Cell::Intersect(Cell* c2, float &i_x, float &i_y){
    float p0_x=p1->y;
    float p0_y=p1->z;
    float p1_x=p2->y;
    float p1_y=p2->z;

    float p2_x=c2->p1->y;
    float p2_y=c2->p1->z;
    float p3_x=c2->p2->y;
    float p3_y=c2->p2->z;


    if(p0_x==p2_x && p0_y==p2_y)return 0;
    if(p1_x==p3_x && p1_y==p3_y)return 0;

    float s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    float s, t;
    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);



    if (s > 0 && s < 1 && t > 0 && t < 1)
    {
        // Collision detected
        i_x = p0_x + (t * s1_x);
        i_y = p0_y + (t * s1_y);
        return 1;
    }

    return 0; // No collision
}

void Cell::PrintLineRepresentation(){

    // x=a+bz

    double b=(p1->y-p2->y)/(p1->z-p2->z);
    double a=-1*b*p1->z+p1->y;
    a+=115*b;
    printf("Cell: (%.3f,%.3f)~(%.3f,%.3f) \t",p1->y,p1->z,p2->y,p2->z);
    printf("Line representation:  %f   %f \n",a, b);

}

void Cell::Print()
{
    printf("Cell:  Layer %i ,%i ,   (%.3f,%.3f)~(%.3f,%.3f),L=%.3f,  State  %f   N Connect Upper: %i \n",p1->layer,p2->layer,p1->y,p1->z,p2->y,p2->z,length,state,v_connected_upper.size());

    //    if(state==11){
    //        cout<<"...........Same Start.........."<<endl;
    //        for(int i=0;i<v_connected_samestart.size();i++){
    //            cout<<"p1 ";  v_connected_samestart.at(i)->p1->Print();
    //            cout<<"p2 ";  v_connected_samestart.at(i)->p2->Print();
    //        }
    //        cout<<"...........Same End.........."<<endl;
    //        for(int i=0;i<v_connected_sameend.size();i++){
    //            cout<<"p1 ";  v_connected_sameend.at(i)->p1->Print();
    //            cout<<"p2 ";  v_connected_sameend.at(i)->p2->Print();
    //        }
    //        cout<<"...........Upper.........."<<endl;
    //            for(int i=0;i<v_connected_lower.size();i++){
    //            cout<<"p1 ";  v_connected_lower.at(i)->p1->Print();
    //            cout<<"p2 ";  v_connected_lower.at(i)->p2->Print();
    //        }
    //        cout<<"...........Lower.........."<<endl;
    //            for(int i=0;i<v_connected_upper.size();i++){
    //                cout<<"p1 ";  v_connected_upper.at(i)->p1->Print();
    //                cout<<"p2 ";  v_connected_upper.at(i)->p2->Print();
    //            }
    //    }

}

void Cell::Draw()
{

    //    TEllipse* point=new TEllipse(hit->x,hit->z,0.5,0.5);

    //    int color=TColor::GetColor(log10(hit->Amp)/log10(5000),float(0.),float(0.));
    //    point->SetFillColor(color);
    //    point->Draw();

    if(state<=0)return;
    double x1,y1,x2,y2;
    x1=p1->y;
    y1=p1->z;
    x2=p2->y;
    y2=p2->z;
    TLine* line=new TLine(x1,y1,x2,y2);
    line->SetLineWidth(state);
    line->SetLineColor(kGreen+2);
    line->Draw();
}

void Cell::ClearConnectedCell(){

    v_connected_lower.clear();
    v_connected_upper.clear();
    v_connected_samestart.clear();
    v_connected_sameend.clear();
    return;

}

void Cell::AddConnectedCell()
{
    ClearConnectedCell();

    for(int i=0;i<p1->v_connected_upper.size();i++){
        if(p1->v_connected_upper.at(i) != this)v_connected_samestart.push_back(p1->v_connected_upper.at(i));
    }

    for(int i=0;i<p2->v_connected_lower.size();i++){
        if(p2->v_connected_lower.at(i) != this)v_connected_sameend.push_back(p2->v_connected_lower.at(i));
    }

    //    double min_angle_lower=999;
    for(int i=0;i<p1->v_connected_lower.size();i++){
        Cell* clower=p1->v_connected_lower.at(i);
        if(clower->state<0)continue;
        //        if(this->GetAbsoluteAngle(clower)<min_angle_lower){
        //        if(IsMatched(clower)) { //    // what is this ?
        //               if(v_connected_lower.size()==0)v_connected_lower.push_back(clower);
        //               else v_connected_lower[0]=clower;
        v_connected_lower.push_back(clower);
        //                min_angle_lower=this->GetAbsoluteAngle(clower);
        //            }
        //        }
    }

    //    double min_angle_upper=999;
    for(int i=0;i<p2->v_connected_upper.size();i++){
        Cell* cupper=p2->v_connected_upper.at(i);
        if(cupper->state<0)continue;
        //        if(this->GetAbsoluteAngle(cupper)<min_angle_upper){
        //        if(IsMatched(cupper)){
        //                if(v_connected_upper.size()==0)v_connected_upper.push_back(cupper);
        //                else v_connected_upper[0]=cupper;
        v_connected_upper.push_back(cupper);
        //                min_angle_upper=this->GetAbsoluteAngle(cupper);
        //            }
        //        }
    }

    return;
}

void Cell::Evolve(){
    if(evolve_next_generation!=0){
        state+=evolve_next_generation;
        evolve_next_generation=0;
    }

}




double Cell::GetTransformFunction(Cell *cell2,double alpha,double gamma)
{
    if(this==cell2)return 0;


    if(this->p1 == cell2->p1  || this->p2 == cell2->p2){  // begin or end in the same point
        return -1*alpha;
    }


    if( (this->p1 != cell2->p2) && (this->p2 != cell2->p1)){ // no common end point
        //            double cell1_mid_x=(p1->x+p2->x)/2;
        //            double cell1_mid_z=(p1->z+p2->z)/2;
        //            double cell2_mid_x=(cell2->p1->x+cell2->p2->x)/2;
        //            double cell2_mid_z=(cell2->p1->z+cell2->p2->z)/2;
        //            double length_between_mid_point=sqrt(pow(cell1_mid_x-cell2_mid_x,2)+pow(cell1_mid_z-cell2_mid_z,2));

        //            if(length_between_mid_point>10)return 0;
        //            double cos_theta_1=fabs(((cell2_mid_x-cell1_mid_x)*(cell2->p2->x-cell2->p1->x)+(cell2_mid_z-cell1_mid_z)*(cell2->p2->z-cell2->p1->z))/(length_between_mid_point*cell2->length));
        //            double cos_theta_2=fabs(((p2->x-p1->x)*(cell2_mid_x-cell1_mid_x)+(p2->z-p1->z)*(cell2_mid_z-cell1_mid_z))/(length_between_mid_point*length));
        //            return 1*pow(cos_theta_1*cos_theta_2,1);


        //        float dummy_x;
        //        float dummy_y;
        //        if(Intersect(cell2,dummy_x,dummy_y)){  // intersect
        //            return -0.5*alpha;
        //        }


        return 0;

    }




    //        double distance=sqrt(pow(p2->x-p1->x+cell2->p2->x-cell2->p1->x,2)+pow(p2->z-p1->z+cell2->p2->z-cell2->p1->z,2));
    double distance=(length+cell2->length);//;/3;


    double cos_theta;  // This is the most correct way,  negative means very large angle
    cos_theta=((p2->y-p1->y)*(cell2->p2->y-cell2->p1->y)+(p2->z-p1->z)*(cell2->p2->z-cell2->p1->z))/(length*cell2->length);

    //    double cos_theta=fabs(cos(GetAbsoluteAngle(cell2)));

    //    double range[2];
    //    GetAngleRange(cell2,range);

    //    double cos_theta_1=cos(range[0]);
    //    double cos_theta_2=cos(range[1]);

    //    //        cout<<"~~~~~~~~~~~~~Angle Range~~~~~~~~~~~~~"<<endl;
    //    //        Print();
    //    //        cell2->Print();
    //    //        cout<<"Angle: "<<range[0]<<", "<<range[1]<<endl;
    //    //        cout<<"cos theta : "<<cos_theta<<endl;
    //    //        cout<<"cos theta 1: "<<cos_theta_1<<endl;
    //    //        cout<<"cos theta 2: "<<cos_theta_2<<endl;

    //    if(range[0]*range[1]<0){
    //        cos_theta=1;
    //    }else if(cos_theta_1>cos_theta_2){
    //        cos_theta=cos_theta_1;
    //    }else cos_theta=cos_theta_2;


    //    cout<<" Final : "<<cos_theta<<endl;

    if(cos_theta<0.9){
        return -1*alpha;
    }



    //    cout<<"Cos Theta: " <<cos_theta<<",  distance: " <<distance<<",   power 54; " <<pow(cos_theta,54)<<endl;


    return gamma*pow(cos_theta,20)/distance;
}



void Cell::InValidate(){
    //    cout<<"Invalidate:  "<<endl;
    //    Print();
    valid=0;
    state=-99;
    return;
}

void Cell::InValidateOtherCells(){


    for(int i=0;i<v_connected_lower.size();i++){
        v_connected_lower.at(i)->InValidate();
    }
    for(int i=0;i<v_connected_upper.size();i++){
        v_connected_upper.at(i)->InValidate();
    }
    //    for(int i=0;i<v_connected_samestart.size();i++){
    //        v_connected_samestart.at(i)->InValidate();
    //    }
    //    for(int i=0;i<v_connected_sameend.size();i++){
    //        v_connected_sameend.at(i)->InValidate();
    //    }


    return;
}

void Cell::SortConnectedCells()
{
    std::sort(v_connected_lower.begin(),v_connected_lower.end(),compCellByCosTheta(this));
    std::sort(v_connected_upper.begin(),v_connected_upper.end(),compCellByCosTheta(this));
    return;
}

int Cell::Merge(Cell *cell2)
{

//    if(cell2->state<0)return 0;
//    cout<<"Merged From: "<<endl;
//    Print();
//    cell2->Print();



//    if(p2==cell2->p1){

//        cout<<"p2==cell2->p1 : "<<endl;

//        for(int i=0;i<v_connected_sameend.size();i++){
//            v_connected_sameend.at(i)->InValidate();
//        }
//        for(int i=0;i<v_connected_upper.size();i++){
//            if(v_connected_upper.at(i) != cell2)v_connected_upper.at(i)->InValidate();
//        }

//        //        p2=cell2->p2;

//        //        v_connected_upper = cell2->v_connected_upper;
//        //        v_connected_sameend = cell2->v_connected_sameend;

//    }else if(p1==cell2->p2){


//        cout<<"p1==cell2->p2 : "<<endl;

//        for(int i=0;i<v_connected_samestart.size();i++){
//            v_connected_samestart.at(i)->InValidate();
//        }
//        for(int i=0;i<v_connected_lower.size();i++){
//            if(v_connected_lower.at(i) != cell2)v_connected_lower.at(i)->InValidate();
//        }

//        //        v_connected_lower = cell2->v_connected_lower;
//        //        v_connected_samestart = cell2->v_connected_samestart;

//        //        p1=cell2->p1;
//    }else{
//        cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
//        cout<<"~~ERROR Merging Not Supported~~~~~~~~"<<endl;
//        cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
//        return 0;
//    }
//    return 0;

//    valid = valid || cell2->valid;

//    HasMatchedCells = HasMatchedCells || cell2->HasMatchedCells;


//    theta=TMath::ATan2((p2->y-p1->y),(p2->z-p1->z));
//    if(p1->type !=p2->type){
//        std::cout<<"~~~~~~Merge~~~~~~~~Error ,  Hit in different projection "<<std::endl;
//        cout<<"p1: " <<p1->layer<<", "<<p2->layer<<endl;
//    }
//    theta_range[0]=TMath::ATan2(((p2->y-0.5)-(p1->y+0.5)),(p2->z-p1->z));
//    theta_range[1]=TMath::ATan2(((p2->y+0.5)-(p1->y-0.5)),(p2->z-p1->z));


//    length=sqrt(pow(p1->y-p2->y,2)+pow(p1->z-p2->z,2));

//    assigned=0;

//    state=1;

//    length=sqrt(pow(p1->y-p2->y,2)+pow(p1->z-p2->z,2));


//    for(int i=0;i<cell2->v_hits.size();i++){
//        if(v_hits[i]!=p1 && v_hits[i]!=p2){
//            v_hits.push_back(v_hits[i]);
//        }
//    }

//    cout<<"Merged To: "<<endl;
//    Print();

//    return 1;

}

double Cell::CosTheta(Cell *cell2)
{
    return ((p2->y-p1->y)*(cell2->p2->y-cell2->p1->y)+(p2->z-p1->z)*(cell2->p2->z-cell2->p1->z))/(length*cell2->length);
}
