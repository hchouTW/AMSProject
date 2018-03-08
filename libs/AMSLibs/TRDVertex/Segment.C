#include "Segment.h"


void Segment::Init()
{
    cells.reserve(20);
    hits.reserve(20);

    x1=0;
    y1=0;
    z1=0;
    x2=0;
    y2=0;
    z2=0;
}

Segment::Segment(){
    Init();
}

Segment::Segment(Cell* start, double theta_threshold)   // This is only for cells constructed with cellular automatum algorithm.
{
    Init();
    Add(start);
    Cell* ptr_current=start;

    //    while(ptr_current->v_connected_upper.size()){
////        ptr_current->Print();
//        Cell* ptr_next=ptr_current->v_connected_upper.at(0);
//        for(int i=1;i<ptr_current->v_connected_upper.size();i++){
//            Cell* ptr_temp=ptr_current->v_connected_upper.at(i);
//            if(ptr_current->state != ptr_temp->state+1)continue;
//            if(ptr_current->GetAbsoluteAngle(ptr_temp)<ptr_current->GetAbsoluteAngle(ptr_next)){
//                ptr_next=ptr_temp;
//            }
//        }
//        if(!ptr_next)break;
//        Add(ptr_next);
//        ptr_current=ptr_next;
//    }

    while(ptr_current->v_connected_upper.size()){
        double costheta_max=theta_threshold;
        double costheta_max_second_choice=theta_threshold;
        Cell* ptr_next=0;
        Cell* ptr_next_second_choice=0;

//        for(int i=0;i<ptr_current->v_connected_upper.size();i++){
//            Cell* ptr_temp=ptr_current->v_connected_upper.at(i);
//            if(ptr_current->state != ptr_temp->state+1)continue;
//            double cos_theta=((ptr_current->p2->y-ptr_current->p1->y)*(ptr_temp->p2->y-ptr_temp->p1->y)+(ptr_current->p2->z-ptr_current->p1->z)*(ptr_temp->p2->z-ptr_temp->p1->z))/(ptr_current->length*ptr_temp->length);
//            if(cos_theta>costheta_max){
//                costheta_max=cos_theta;
//                ptr_next=ptr_temp;
//            }
//        }
//        if(!ptr_next)break;
//        Add(ptr_next);
//        ptr_current=ptr_next;
//    }


        for(int i=0;i<ptr_current->v_connected_upper.size();i++){
            Cell* ptr_temp=ptr_current->v_connected_upper.at(i);
            double cos_theta=((ptr_current->p2->y-ptr_current->p1->y)*(ptr_temp->p2->y-ptr_temp->p1->y)+(ptr_current->p2->z-ptr_current->p1->z)*(ptr_temp->p2->z-ptr_temp->p1->z))/(ptr_current->length*ptr_temp->length);
            if(ptr_current->state == ptr_temp->state+1){
                if(cos_theta>costheta_max){
                    costheta_max=cos_theta;
                    ptr_next=ptr_temp;
                }
            }else{
                if(cos_theta>costheta_max_second_choice){
                    costheta_max_second_choice=cos_theta;
                    ptr_next_second_choice=ptr_temp;
                }
            }
        }

        if(!ptr_next)ptr_next=ptr_next_second_choice;
        if(!ptr_next)break;
        Add(ptr_next);
        ptr_current=ptr_next;
    }




//    Finish();
}



int Segment::Add(Cell *_cell){
    if(_cell->state<0)return 0;
    int isfirst=cells.size()==0;
    if(isfirst){
        cells.push_back(_cell);
        //            _cell->p1->Print();
        //            _cell->p2->Print();
        hits.push_back(_cell->p1);
        hits.push_back(_cell->p2);

        x1=_cell->p1->y;
        z1=_cell->p1->z;
    }else{
        Cell* last=cells.at(cells.size()-1);
        if(last->p2!=_cell->p1){
            return 0;
        }
        if(_cell->state!=last->state-1){
            return 0;
        }
        if(!last->IsMatched(_cell)){
            return 0;
        }
        //            _cell->p2->Print();
        cells.push_back(_cell);
        hits.push_back(_cell->p2);
        x2=_cell->p2->y;
        z2=_cell->p2->z;
    }

    return 1;
}

int Segment::Finish(){
//    std::cout<<"Finish: " <<std::endl;
    for(int i=0;i<cells.size();i++){
        (cells.at(i))->assigned++;
//        cells.at(i)->Print();
    }

}

void Segment::Print(){
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"Segment with "<<GetNHits()<<" hits  and  "<<GetNCells()<<" Cells"<<std::endl;
    std::cout<<"Track Parameter: "<<slope<<",   "<<intercept<<std::endl;
    for(int i=0;i<cells.size();i++){
        std::cout<<"\t";cells.at(i)->Print();
    }
    std::cout<<"------------------------------"<<std::endl;

}



void Segment::FitStraightLine(double& chi2, double& ndof){
    double xmean=0;
    double ymean=0;
    double N=0;
    for(int i=0;i<hits.size();i++){
        TRD2DHit* hit=hits.at(i);
//                    hit->Print();
        ymean+=hit->z;
        xmean+=hit->y;
        N++;
    }
    xmean/=N;
    ymean/=N;
    double s_x=0;
    double s_y=0;
    double s_xy=0;
    double s_yy=0;
    for(int i=0;i<hits.size();i++){
        TRD2DHit* hit=hits.at(i);
        s_y=hit->z-ymean;
        s_x=hit->y-xmean;
        s_xy+=s_x*s_y;
        s_yy+=s_y*s_y;
    }
    if(s_xy!=0)slope=s_yy/s_xy;    // Fit on "x" residule, then invert the slope to get normal X-Y convention
    else slope=1e10;
    intercept=ymean-slope*xmean;



    // Simple Chi2
    chi2=0;
    ndof=hits.size();
    for(int i=0;i<hits.size();i++){
        TRD2DHit* hit=hits.at(i);
        chi2+=pow((hit->y-(hit->z-intercept)/slope)/0.3,2);
    }
    ndof-=2;
    return;
}

void Segment::Draw(){
    double plot_y1=-200;
    double plot_y2=200;
    double plot_x1=(plot_y1-intercept)/slope;
    double plot_x2=(plot_y2-intercept)/slope;
    TLine* line=new TLine(plot_x1,plot_y1,plot_x2,plot_y2);
//    std::cout<<"P1: " <<plot_x1<<", "<<plot_y1<<std::endl;
//    std::cout<<"P2: " <<plot_x2<<", "<<plot_y2<<std::endl;
    line->SetLineWidth(2);
    line->SetLineColor(kRed);
    line->Draw();
}
