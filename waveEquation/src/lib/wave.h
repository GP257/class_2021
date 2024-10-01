#include<float3DReg.h>
#include<float1DReg.h>

using namespace SEP;
class acousticDensity{
	public:
 
		acousticDensity(){;}
   acousticDensity(std::shared_ptr<float3DReg> density,
		   std::shared_ptr<float3DReg> lamda,
		   const int ix, const int iy, const int iz,
		   std::shared_ptr<float1DReg> sval,
		   const float fmax,
		   const float dt,
		   const int islice,
		   std::vector<std::vector<int>> bV,
		   std::vector<std::vector<int>> eV,
		   std::vector<std::vector<int>> bS,
		   std::vector<std::vector<int>> eS){
          basicInit(density,lamda,ix,iy,iz,sval,fmax,dt,islice,bV,eV,bS,eS);

         }
     void basicInit(std::shared_ptr<float3DReg> density,
		   std::shared_ptr<float3DReg> lamda,
		   const int ix, const int iy, const int iz,
		   std::shared_ptr<float1DReg> sval,
		   const float fmax,
		   const float dt,
		   const int islice,
           std::vector<std::vector<int>> bV,
		   std::vector<std::vector<int>> eV,
		   std::vector<std::vector<int>> bS,
		   std::vector<std::vector<int>> eS);


	
   void prop(const int nt);
   std::shared_ptr<float3DReg> getWavefield(){
	   return _wave;
  }
   virtual void calcStress();
   virtual void calcVel();

  protected:
     int _nts;
    float _dt;
    int _ix,_iy,_iz,_islice;
    std::shared_ptr<float1DReg> _source;
   float _x[5],_y[5],_z[5];
   std::shared_ptr<float3DReg> _Vx,_Vy,_Vz,_prev,_cur,_next,_lamda,_invD,_wave;
   std::vector<int> _b1V,_e1V,_b1S,_e1S;
   std::vector<int> _b2V,_e2V,_b2S,_e2S;
   std::vector<int> _b3V,_e3V,_b3S,_e3S,_ns;
};
class acousticDensityISPC: public acousticDensity{
	public:


   acousticDensityISPC(std::shared_ptr<float3DReg> density,
		   std::shared_ptr<float3DReg> lamda,
		   const int ix, const int iy, const int iz,
		   std::shared_ptr<float1DReg> sval,
		   const float fmax,
		   const float dt,
		   const int islice,
		   std::vector<std::vector<int>> bV,
		   std::vector<std::vector<int>> eV,
		   std::vector<std::vector<int>> bS,
		   std::vector<std::vector<int>> eS){
         basicInit(density,lamda,ix,iy,iz,sval,fmax,dt,islice,bV,eV,bS,eS);

         }


   virtual void calcStress() override;
   virtual void calcVel() override;

};

class acousticDensityISPC_P: public acousticDensity{
	public:


   acousticDensityISPC_P(std::shared_ptr<float3DReg> density,
		   std::shared_ptr<float3DReg> lamda,
		   const int ix, const int iy, const int iz,
		   std::shared_ptr<float1DReg> sval,
		   const float fmax,
		   const float dt,
		   const int islice,
		   std::vector<std::vector<int>> bV,
		   std::vector<std::vector<int>> eV,
		   std::vector<std::vector<int>> bS,
		   std::vector<std::vector<int>> eS){
         basicInit(density,lamda,ix,iy,iz,sval,fmax,dt,islice,bV,eV,bS,eS);

         }


   virtual void calcStress() override;
   virtual void calcVel() override;

};
