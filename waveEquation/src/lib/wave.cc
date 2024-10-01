#include<wave.h>
#include "kernel.h"
#include <omp.h>

void acousticDensity::basicInit(std::shared_ptr<float3DReg> density,
		   std::shared_ptr<float3DReg> lamda,
		   const int ix, const int iy, const int iz,
		   std::shared_ptr<float1DReg> sval, const float fmax,
		   const float dt,
		   const int islice,
		   std::vector<std::vector<int>> bV,
		   std::vector<std::vector<int>> eV,
		   std::vector<std::vector<int>> bS,
		   std::vector<std::vector<int>> eS){

    float minV=(*lamda->_mat)[0][0][0]/(*density->_mat)[0][0][0];
    float maxV=minV;

    std::vector<axis> axes=lamda->getHyper()->getAxes();
    float dmax=std::max(std::max(axes[0].d,axes[1].d),axes[2].d);
    float dmin=std::min(std::min(axes[0].d,axes[1].d),axes[2].d);

    std::vector<int> ns=lamda->getHyper()->getNs();

    _ix=ix;
    _iy=iy;
    _iz=iz;
    _dt=dt;
    _source=sval;
    _b3V=bV[0];
    _b2V=bV[1];
    _b1V=bV[2];
    _b3S=bS[0];
    _b2S=bS[1];
    _b1S=bS[2];
    _e3V=eV[0];
    _e2V=eV[1];
    _e1V=eV[2];
    _e3S=eS[0];
    _e2S=eS[1];
    _e1S=eS[2];
    _islice=islice;
    
     _ns=lamda->getHyper()->getNs();

    _nts=sval->getHyper()->getAxes()[0].n;
    this->_dt=dt;
    this->_invD=lamda->clone();
    this->_lamda=lamda;
    this->_Vx=lamda->clone();
    this->_Vy=(lamda->clone());
    this->_Vz=(lamda->clone());
    this->_prev=(lamda->clone());
    this->_cur=(lamda->clone());
    this->_next=(lamda->clone());

    for(int i3=0; i3 < ns[2]; i3++){
	    for(int i2=0; i2 < ns[1]; i2++){
		    for(int i1=0; i1 < ns[0]; i1++){
			    (*_invD->_mat)[i3][i2][i1]=(*density->_mat)[i3][i2][i1];
			    maxV=std::max(maxV,(*lamda->_mat)[i3][i2][i1]/(*density->_mat)[i3][i2][i1]);
			    minV=std::min(minV,(*lamda->_mat)[i3][i2][i1]/(*density->_mat)[i3][i2][i1]);
			}
	    }
    }
    if(maxV*dt /dmin >.5) throw SEPException("Unstable");
    std::cerr<<"cbeck "<<minV<<" "<<fmax<<" "<<dmax<<std::endl;
    if(minV/fmax/dmax<2.8) throw SEPException("Dispeserive");

    float co[45]={-1.18679470e-04, 1.76565988e-03,-1.38427734e-02, 8.97216797e-02,-1.21124268e+00};
    for(int i=0;i < 5; i++){
	    _x[i]=co[i]/axes[0].d;
	    _y[i]=co[i]/axes[0].d;
	    _z[i]=co[i]/axes[0].d;
    }
    _prev->zero();
    _cur->zero();
    _next->zero();
}



	
void acousticDensity::prop(const int nt){
	std::vector<axis> axes=_lamda->getHyper()->getAxes();
        axes[2].n=nt;
	std::shared_ptr<float3DReg> x(new float3DReg(axes[0],axes[1],axes[2]));
	_wave=x;
	for(int it=0;it < nt; it++){
		std::cerr<<"looping "<<it<<std::endl;
           calcVel();
	   calcStress();
     if(it<_nts)
	   (*_next->_mat)[_ix][_iy][_iz]+=(*_source->_mat)[it];


	   for(int i2=0; i2 < axes[1].n; i2++){
		   for(int i1=0; i1 < axes[0].n;i1++){
	   (*_wave->_mat)[it][i2][i1]=(*_next->_mat)[_islice][i2][i1];
	   }
		   }
           
	    std::shared_ptr<float3DReg> t=_prev;
	    _prev=_cur;
	    _cur=_next;
       }

	   }



void acousticDensity::calcStress(){

	
for(int ib=0;ib< _b3S.size(); ib++){
  for(int i3=_b3S[ib]; i3 < _e3S[ib];i3++){
    for(int i2=_b2S[ib]; i2 <_e2S[ib];i2++){
      for(int i1=_b1S[ib]; i1<_e1S[ib];i1++){
	 float x=(
          _x[4]*((*_Vx->_mat)[i3][i2][i1-0]-(*_Vx->_mat)[i3][i2][i1+1])+
          _x[3]*((*_Vx->_mat)[i3][i2][i1-1]-(*_Vx->_mat)[i3][i2][i1+2])+
          _x[2]*((*_Vx->_mat)[i3][i2][i1-2]-(*_Vx->_mat)[i3][i2][i1+3])+
          _x[1]*((*_Vx->_mat)[i3][i2][i1-3]-(*_Vx->_mat)[i3][i2][i1+4])+
          _x[0]*((*_Vx->_mat)[i3][i2][i1-4]-(*_Vx->_mat)[i3][i2][i1+5]));
	 float y=
          _y[4]*((*_Vy->_mat)[i3][i2-0][i1]-(*_Vy->_mat)[i3][i2+1][i1])+
          _y[3]*((*_Vy->_mat)[i3][i2-1][i1]-(*_Vy->_mat)[i3][i2+2][i1])+
          _y[2]*((*_Vy->_mat)[i3][i2-2][i1]-(*_Vy->_mat)[i3][i2+3][i1])+
          _y[1]*((*_Vy->_mat)[i3][i2-3][i1]-(*_Vy->_mat)[i3][i2+4][i1])+
          _y[0]*((*_Vy->_mat)[i3][i2-4][i1]-(*_Vy->_mat)[i3][i2+5][i1]);
	  float z=
          _z[4]*((*_Vz->_mat)[i3][i2][i1]-(*_Vz->_mat)[i3+1][i2][i1])+
          _z[3]*((*_Vz->_mat)[i3-1][i2][i1]-(*_Vz->_mat)[i3+2][i2][i1])+
          _z[2]*((*_Vz->_mat)[i3-2][i2][i1]-(*_Vz->_mat)[i3+3][i2][i1])+
          _z[1]*((*_Vz->_mat)[i3-3][i2][i1]-(*_Vz->_mat)[i3+4][i2][i1])+
          _z[0]*((*_Vz->_mat)[i3-4][i2][i1]-(*_Vz->_mat)[i3+5][i2][i1]);
          double right= (*_invD->_mat)[i3][i2][i1]*(x+y+z);
	  (*_next->_mat)[i3][i2][i1]=2*(*_cur->_mat)[i3][i2][i1]-(*_prev->_mat)[i3][i2][i1]+_dt*_dt*right;

      }}}}
}
void acousticDensity::calcVel(){

for(int ib=0;ib< _b3V.size(); ib++){
  for(int i3=_b3V[ib]; i3 < _e3V[ib];i3++){
    for(int i2=_b2V[ib]; i2 <_e2V[ib];i2++){
      for(int i1=_b1V[ib]; i1<_e1V[ib];i1++){
        (*_Vx->_mat)[i3][i2][i1]=(*_lamda->_mat)[i3][i2][i1]*(
          _x[4]*((*_cur->_mat)[i3][i2][i1-1]-(*_cur->_mat)[i3][i2][i1])+
          _x[3]*((*_cur->_mat)[i3][i2][i1-2]-(*_cur->_mat)[i3][i2][i1+1])+
          _x[2]*((*_cur->_mat)[i3][i2][i1-3]-(*_cur->_mat)[i3][i2][i1+2])+
          _x[1]*((*_cur->_mat)[i3][i2][i1-4]-(*_cur->_mat)[i3][i2][i1+3])+
          _x[0]*((*_cur->_mat)[i3][i2][i1-5]-(*_cur->_mat)[i3][i2][i1+4]));
         (*_Vy->_mat)[i3][i2][i1]=(*_lamda->_mat)[i3][i2][i1]*(
          _y[4]*((*_cur->_mat)[i3][i2-1][i1]-(*_cur->_mat)[i3][i2][i1])+
          _y[3]*((*_cur->_mat)[i3][i2-2][i1]-(*_cur->_mat)[i3][i2+1][i1])+
          _y[2]*((*_cur->_mat)[i3][i2-3][i1]-(*_cur->_mat)[i3][i2+2][i1])+
          _y[1]*((*_cur->_mat)[i3][i2-4][i1]-(*_cur->_mat)[i3][i2+3][i1])+
          _y[0]*((*_cur->_mat)[i3][i2-5][i1]-(*_cur->_mat)[i3][i2+4][i1]));
         (*_Vz->_mat)[i3][i2][i1]=(*_lamda->_mat)[i3][i2][i1]*(
          _z[4]*((*_cur->_mat)[i3-1][i2][i1]-(*_cur->_mat)[i3][i2][i1])+
          _z[3]*((*_cur->_mat)[i3-2][i2][i1]-(*_cur->_mat)[i3+1][i2][i1])+
          _z[2]*((*_cur->_mat)[i3-3][i2][i1]-(*_cur->_mat)[i3+2][i2][i1])+
          _z[1]*((*_cur->_mat)[i3-4][i2][i1]-(*_cur->_mat)[i3+3][i2][i1])+
          _z[0]*((*_cur->_mat)[i3-5][i2][i1]-(*_cur->_mat)[i3+4][i2][i1]));
      }}}}
}


void acousticDensityISPC::calcStress(){

   stressISPC(_b3S.size(),&_ns[0],&_b3V[0],&_e3V[0],&_b2V[0],&_e2V[0],&_b1V[0],
   &_e1V[0],
   &_x[0],&_y[0],&_z[0],
   _dt*_dt,
   _Vx->getVals(),_Vy->getVals(),_Vz->getVals(),
   _invD->getVals(),_prev->getVals(),_cur->getVals(),
   _next->getVals());


}
void acousticDensityISPC_P::calcVel() {

#pragma omp parallel for 
  for(int ib=0; ib < _b3S.size();ib++){
   velISPC_B(&_ns[0],_b3V[ib],_e3V[ib],_b2V[ib],_e2V[ib],_b1V[ib],
   _e1V[ib],
   &_x[0],&_y[0],&_z[0],
   _cur->getVals(),_lamda->getVals(),
   _Vx->getVals(),_Vy->getVals(),_Vz->getVals());
  }

}


void acousticDensityISPC_P::calcStress(){


#pragma omp parallel for 
  for(int ib=0; ib < _b3S.size();ib++){
   stressISPC_B(&_ns[0],_b3V[ib],_e3V[ib],_b2V[ib],_e2V[ib],_b1V[ib],
   _e1V[ib],
   &_x[0],&_y[0],&_z[0],
   _dt*_dt,
   _Vx->getVals(),_Vy->getVals(),_Vz->getVals(),
   _invD->getVals(),_prev->getVals(),_cur->getVals(),
   _next->getVals());
  }


}
void acousticDensityISPC::calcVel() {

   velISPC(_b3V.size(),&_ns[0],&_b3V[0],&_e3V[0],&_b2V[0],&_e2V[0],&_b1V[0],
   &_e1V[0],
   &_x[0],&_y[0],&_z[0],
   _cur->getVals(),_lamda->getVals(),
   _Vx->getVals(),_Vy->getVals(),_Vz->getVals());

}



