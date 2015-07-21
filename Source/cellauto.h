/*cellauto.h contains classes and operations pertaining to the construction,
and manipulation of cellular automata.

vector,math.h,time.h, and GL/glut headers must be included before including
this header file*/

#ifndef __CELLAUTO_H_INCLUDED__   // if cellauto.h hasn't been included yet...
#define __CELLAUTO_H_INCLUDED__   //   #define this so the compiler knows it has been included
#include <math.h>
#include <vector>
using namespace std;

class cellauto{
    int rank,size,interrad,resetoption,currule;
    vector<int> rule;
    vector<short> cell,walsh; //unformatted cellauto data
    vector<vector<short> > cellhist,walshhist; //record of past cellauto configurations
    
    public:
    int getsize(){return size;}
    void setinterrad(int ir){interrad=ir;} //set interaction radius.  Default is 1
    int getinterrad(){return interrad;} //return value of interaction radius
    void rulechange(short direction){
        currule=currule+direction;
        if(currule==rule.size())randrule();
        else if(currule<1)currule=0;
        cout<<"\ncurrent rule: "<<rule[currule];
        reset(resetoption);
    }
     
    void reset(int option){
        resetoption=option;
        cell.clear();
        cellhist.clear();
        walsh.clear();
        walshhist.clear();
        int i,halfway;
        halfway=(((1-pow(size,rank+1))/(1-size))-1)/2; //index of point that sits in the middle of the hypercube of given size and rank
        switch(option){
            case 1: for(i=0;i<pow(size,rank);i++){//precalculate 2 frames using delta function at space center
                        if(i==halfway)cell.push_back(1); //set central element = 1
                        else cell.push_back(0); //set all other elements = 0
                    }
                    cellhist.push_back(cell);
                    cellhist.push_back(cell);
                    break;
            case 2: for(i=0;i<pow(size,rank);i++)cell.push_back(rand()%2);//precalculate 2 frames using random values
                    cellhist.push_back(cell);
                    cell.clear();
                    for(i=0;i<pow(size,rank);i++)cell.push_back(rand()%2);
                    cellhist.push_back(cell);
                    break;
        }
    }
    void randrule(){ //set new random rule
        rule.push_back(rand()%(int)pow(2,(pow(((interrad*2)+1),rank)+1)));
        reset(resetoption);
    }
    
    void init(int r,int s,int ir,int ruleseed){
        srand(time(NULL));
        currule=0;
        rank=r; //define rank of CA
        size=pow(2,floor(log(s)/log(2))); //define size of CA by rounding down to nearest power of 2
        interrad=ir; //specified interaction radius
        rule.push_back(ruleseed%(int)pow(2,(pow(((interrad*2)+1),rank)+1))); //set random rule
        cout<<"\ncurrent rule: "<<rule[currule];
        resetoption=1; //default reset type (1=delta, 2=random)
        reset(resetoption);
    }
    
    cellauto(int r,int s,int ir,int ruleseed){
        init(r,s,ir,ruleseed);
    }
    cellauto(){}
    
    void cellrecord(){
        if(cellhist.size()==size)cellhist.erase(cellhist.begin()+0);
        int i,j,neighborloc,totalpos,curstate;
        vector<short> temprelpos,relpos; //relative position in [rank] dimensional space
        for(i=0;i<rank;i++)temprelpos.push_back(-interrad); //start at neg. interaction radius
        for(i=0;i<pow(size,rank);i++){
            relpos=temprelpos;
            totalpos=0; //0->N index representing current neighbor
            curstate=0; //determined by a sum over values of neighbors within interaction region
            while(relpos[rank-1]<=interrad){ //halt when pos. in last dimension > interaction radius
                neighborloc=i;
                for(j=0;j<rank;j++)neighborloc+=relpos[j]*pow(size,j); //convert multi-dimensional displacement to 1D offset
                if(neighborloc<0 || neighborloc>=pow(size,rank))neighborloc=(int)(neighborloc+pow(size,rank))%(int)pow(size,rank);
                curstate+=cellhist.back()[neighborloc]*pow(2,totalpos); //factor in surrounding neighbors of cell
                for(j=0;j<rank;j++){
                    if(relpos[j]==interrad){
                        if(j==rank-1)relpos[j]++;
                        else{
                            relpos[j]=-interrad;
                            relpos[j+1]++;
                        }
                    }
                    else{
                        relpos[j]++;
                        break;
                    }
                }
            }
            curstate+=cellhist[cellhist.size()-2][i]*pow(2,totalpos); //factor in 2nd order past state of current cell
            cell[i]=(rule[currule] >> curstate) & 1; //set new cell state based on curstate
        }
        cellhist.push_back(cell);
    };
    
    void timewalshrecord(){ //perform walsh-hadamard transform and store result in both walsh and walshhist
        short sign;
        vector<short> tempwalsh;
        int i,j,k;
        int halfway=(((1-pow(size,rank+1))/(1-size))-1)/2; //index of point that sits in the middle of the hypercube of given size and rank
        for(i=halfway-ceil(size/2);i<halfway+floor(size/2);i++)tempwalsh.push_back(cell[i]); //extract line of data passing through central point of hypercube
        
        walsh=cell;
        if(walshhist.size()==size)walshhist.erase(walshhist.begin());
        for(i=0;i<pow(size,rank-1);i+=size){
            for(j=size/2;j>=1;j/=2){
                sign=-1;
                for(k=0;k<size;k++){
                    if(k%j==0)sign*=-1;
                    walsh[k+i]=(sign*tempwalsh[k])+tempwalsh[k+(sign*j)];
                }
                for(k=0;k<size;k++)tempwalsh[k]=walsh[k+i];
            }
        }
        walshhist.push_back(walsh);
    }
    
    void timedisp(){
        int i,j;
        for(j=0;j<cellhist.size();j++)
            for(i=0;i<size;i++){
                glColor3f(cellhist[j][i],cellhist[j][i],cellhist[j][i]);
                glBegin(GL_QUADS);
                    glVertex2f((float)(i%size)/size,(float)j/size);
                    glVertex2f((float)((i%size)+1)/size,(float)j/size);
                    glVertex2f((float)((i%size)+1)/size,(float)(j+1)/size);
                    glVertex2f((float)(i%size)/size,(float)(j+1)/size);
                glEnd();
            }
    }
    
    void spacedisp(){
        int i,visrank;
        if(rank==1)visrank=1;
        else visrank=2;
        
        for(i=0;i<pow(size,visrank);i++){
            glColor3f(cell[i],cell[i],cell[i]);
            glBegin(GL_QUADS);
                glVertex2f((float)(i%size)/size,(float)floor(i/size)/size);
                glVertex2f((float)((i%size)+1)/size,(float)floor(i/size)/size);
                glVertex2f((float)((i%size)+1)/size,(float)(floor(i/size)+1)/size);
                glVertex2f((float)(i%size)/size,(float)(floor(i/size)+1)/size);
            glEnd();
        }
    }
    
    void groupsumdisp(short xrad,short yrad,short sqrad){
        int i,j,k,l,visrank;
        float xtot,ytot,sqtot;
        if(rank==1){
        for(i=0;i<cellhist.size();i++)
            for(j=0;j<size;j++){
                xtot=0;ytot=0;sqtot=0;
                for(k=-xrad;k<=xrad;k++)xtot+=cellhist[i][(k+j+size)%size];
                for(k=-yrad;k<=yrad;k++)ytot+=cellhist[(k+i+cellhist.size())%cellhist.size()][j];
                for(k=-sqrad;k<=sqrad;k++)
                    for(l=-sqrad;l<=sqrad;l++)sqtot+=cellhist[(i+l+cellhist.size())%cellhist.size()][(j+k+size)%size];
            glColor3f((xrad>0)*xtot/(2*xrad+1),(yrad>0)*ytot/(2*yrad+1),(sqrad>0)*sqtot/pow(2*sqrad+1,2));
            glBegin(GL_QUADS);
                glVertex2f((float)(j%size)/size,(float)(i%cellhist.size())/size);
                glVertex2f((float)((j%size)+1)/size,(float)(i%cellhist.size())/size);
                glVertex2f((float)((j%size)+1)/size,(float)((i%cellhist.size())+1)/size);
                glVertex2f((float)(j%size)/size,(float)((i%cellhist.size())+1)/size);
            glEnd();
            }
        }
        else if(rank>1){
        for(i=0;i<size;i++)
            for(j=0;j<size;j++){
                xtot=0;ytot=0;sqtot=0;
                for(k=-xrad;k<=xrad;k++)xtot+=cell[(k+j+size)%size+size*i];
                for(k=-yrad;k<=yrad;k++)ytot+=cell[j+(int)(size*(k+i)+pow(size,2))%(int)pow(size,2)];
                for(k=-sqrad;k<=sqrad;k++)
                    for(l=-sqrad;l<=sqrad;l++)sqtot+=cell[(j+k)%size+(size*(i+l)+(int)pow(size,2))%(int)pow(size,2)];
                glColor3f((xrad>0)*xtot/(2*xrad+1),(yrad>0)*ytot/(2*yrad+1),(sqrad>0)*sqtot/pow(2*sqrad+1,2));
                glBegin(GL_QUADS);
                    glVertex2f((float)(j%size)/size,(float)(i%size)/size);
                    glVertex2f((float)((j%size)+1)/size,(float)(i%size)/size);
                    glVertex2f((float)((j%size)+1)/size,(float)((i%size)+1)/size);
                    glVertex2f((float)(j%size)/size,(float)((i%size)+1)/size);
                glEnd();
            }
        }
    }
    void groupproddisp(short xrad,short yrad,short sqrad){
        int i,j,k,l,visrank;
        float xtot,ytot,sqtot;
        if(rank==1){
        for(i=0;i<cellhist.size();i++)
            for(j=0;j<size;j++){
                xtot=1;ytot=1;sqtot=1;
                for(k=-xrad;k<=xrad;k++)xtot*=cellhist[i][(k+j+size)%size];
                for(k=-yrad;k<=yrad;k++)ytot*=cellhist[(k+i+cellhist.size())%cellhist.size()][j];
                for(k=-sqrad;k<=sqrad;k++)
                    for(l=-sqrad;l<=sqrad;l++)sqtot*=cellhist[(i+l+cellhist.size())%cellhist.size()][(j+k+size)%size];
            glColor3f((xrad>0)*xtot,(yrad>0)*ytot,(sqrad>0)*sqtot);
            glBegin(GL_QUADS);
                glVertex2f((float)(j%size)/size,(float)(i%cellhist.size())/size);
                glVertex2f((float)((j%size)+1)/size,(float)(i%cellhist.size())/size);
                glVertex2f((float)((j%size)+1)/size,(float)((i%cellhist.size())+1)/size);
                glVertex2f((float)(j%size)/size,(float)((i%cellhist.size())+1)/size);
            glEnd();
            }
        }
        else if(rank>1){
        for(i=0;i<size;i++)
            for(j=0;j<size;j++){
                xtot=1;ytot=1;sqtot=1;
                for(k=-xrad;k<=xrad;k++)xtot*=cell[(k+j+size)%size+size*i];
                for(k=-yrad;k<=yrad;k++)ytot*=cell[j+(int)(size*(k+i)+pow(size,2))%(int)pow(size,2)];
                for(k=-sqrad;k<=sqrad;k++)
                    for(l=-sqrad;l<=sqrad;l++)sqtot*=cell[(j+k)%size+(size*(i+l)+(int)pow(size,2))%(int)pow(size,2)];
                glColor3f((xrad>0)*xtot,(yrad>0)*ytot,(sqrad>0)*sqtot);
                glBegin(GL_QUADS);
                    glVertex2f((float)(j%size)/size,(float)(i%size)/size);
                    glVertex2f((float)((j%size)+1)/size,(float)(i%size)/size);
                    glVertex2f((float)((j%size)+1)/size,(float)((i%size)+1)/size);
                    glVertex2f((float)(j%size)/size,(float)((i%size)+1)/size);
                glEnd();
            }
        }
    }
    
    void freqdisp(){
        int i,j;
        for(j=0;j<walshhist.size();j++)
            for(i=0;i<size;i++){
                glColor3f(walshhist[j][i],walshhist[j][i],walshhist[j][i]);
                glBegin(GL_QUADS);
                    glVertex2f((float)(i%size)/size,(float)j/size);
                    glVertex2f((float)((i%size)+1)/size,(float)j/size);
                    glVertex2f((float)((i%size)+1)/size,(float)(j+1)/size);
                    glVertex2f((float)(i%size)/size,(float)(j+1)/size);
                glEnd();
            }
    }
};

#endif 
