#include<wave.h>
int main(int argc, char **argv){
	std::cerr<<"what 1"<<std::endl;
	std::vector<axis> axes;
	axes.push_back(axis(200,0.,.005));
	axes.push_back(axis(200,0.,.005));
	axes.push_back(axis(200,0.,.005));
	std::shared_ptr<hypercube> hyper(new hypercube(axes));
	std::shared_ptr<float3DReg> dens(new float3DReg(hyper)); 
	std::cerr<<"2hat 1"<<std::endl;
	std::shared_ptr<float3DReg> lamda(new float3DReg(hyper)); 
	std::cerr<<"3hat 1"<<std::endl;
	std::shared_ptr<float1DReg> sval(new float1DReg(axes[0]));
	dens->set(1.);
	lamda->set(1.);
	std::cerr<<"4hatx1"<<std::endl;
	sval->zero();
	std::cerr<<"5hat v"<<std::endl;
	(*sval->_mat)[0]=1.;
	std::vector<std::vector<int>> bV,eV,bS,eS;
	std::cerr<<"what 1"<<std::endl;
	std::vector<int> b1,b2,b3,e1,e2,e3;
	b1.push_back(5);
	b2.push_back(5);
	b3.push_back(5);
	e1.push_back(195);
	e2.push_back(195);
	e3.push_back(195);
	bV.push_back(b1);
	bV.push_back(b2);
	bV.push_back(b3);
	eV.push_back(e1);
	eV.push_back(e2);
	eV.push_back(e3);
	bS=bV;
	eS=eV;
	std::cerr<<"AAA 1"<<std::endl;

       acousticDensityISPC x=acousticDensityISPC(dens,lamda,100,100,100,sval,20.,.001,100,bV,eV,bS,eS);
	std::cerr<<"BB 1"<<std::endl;
       
       x.prop(100);

}
