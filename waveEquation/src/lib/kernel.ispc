export void  stressISPC(uniform int nb,uniform int ns[],uniform  int b3[],uniform int e3[],
uniform int b2[], uniform int e2[], uniform int b1[],
uniform int e1[],
uniform float x[], uniform float y[],uniform float z[],
uniform float dt2,
uniform float Vx[], uniform float Vy[], uniform float Vz[],
uniform float invD[], uniform float prev[], uniform float cur[],
uniform float nex[]){

uniform int n1=ns[0],n12=ns[0]*ns[1];

uniform int ib,i2,i3;
for( ib=0;ib < nb; ib++){
  for( i3=b3[ib]; i3 < e3[ib]; i3++){
    for( i2=b2[ib]; i2 < e2[ib]; i2++){
      foreach(i1=b1[ib]... e1[ib]){
        int i=i1+i2*n1+i3*n12;
         nex[i]=2*cur[i]-prev[i]+dt2*invD[i]*(
           x[4]*(Vx[i-0]-Vx[i+1])+
           x[3]*(Vx[i-1]-Vx[i+2])+
           x[2]*(Vx[i-2]-Vx[i+3])+
           x[1]*(Vx[i-3]-Vx[i+4])+
           x[0]*(Vx[i-4]-Vx[i+5])+
           y[4]*(Vy[i-0*n1]-Vy[i+0])+
           y[3]*(Vy[i-1*n1]-Vy[i+2*n1])+
           y[2]*(Vy[i-2*n1]-Vy[i+3*n1])+
           y[1]*(Vy[i-3*n1]-Vy[i+4*n1])+
           y[0]*(Vy[i-4*n1]-Vy[i+5*n1])+
           z[4]*(Vz[i-0*n12]-Vz[i+0])+
           z[3]*(Vz[i-1*n12]-Vz[i+2*n12])+
           z[2]*(Vz[i-2*n12]-Vz[i+3*n12])+
           z[1]*(Vz[i-3*n12]-Vz[i+4*n12])+
           z[0]*(Vz[i-4*n12]-Vz[i+5*n12]));
       }
    }
  }
        
}
}

export void  velISPC(uniform int nb, uniform int ns[],uniform int b3[],uniform int e3[],
uniform int b2[], uniform int e2[], uniform int b1[],
uniform int e1[],
uniform float x[], uniform float y[],uniform float z[],
uniform float cur[], uniform float lamda[],
uniform float Vx[], uniform float Vy[], uniform float Vz[]){

uniform int n1=ns[0],n12=ns[0]*ns[1];
uniform int ib,i2,i3;
for( ib=0;ib < nb; ib++){
  for( i3=b3[ib]; i3 < e3[ib]; i3++){
    for( i2=b2[ib]; i2 < e2[ib]; i2++){
      foreach(i1=b1[ib]... e1[ib]){
        int i=i1+i2*n1+i3*n12;
        Vx[i]=lamda[i]*(
           x[4]*(cur[i-1]-cur[i+0])+
           x[3]*(cur[i-2]-cur[i+1])+
           x[2]*(cur[i-3]-cur[i+2])+
           x[1]*(cur[i-4]-cur[i+3])+
           x[0]*(cur[i-5]-cur[i+4]));

        Vy[i]=lamda[i]*(
           y[4]*(cur[i-1*n1]-cur[i+0])+
           y[3]*(cur[i-2*n1]-cur[i+1*n1])+
           y[2]*(cur[i-3*n1]-cur[i+2*n1])+
           y[1]*(cur[i-4*n1]-cur[i+3*n1])+
           y[0]*(cur[i-5*n1]-cur[i+4*n1]));

        Vz[i]=lamda[i]*(
           z[4]*(cur[i-1*n12]-cur[i+0])+
           z[3]*(cur[i-2*n12]-cur[i+1*n12])+
           z[2]*(cur[i-3*n12]-cur[i+2*n12])+
           z[1]*(cur[i-4*n12]-cur[i+3*n12])+
           z[0]*(cur[i-5*n12]-cur[i+4*n12]));
        
        }
     }
   }
 }
}
export void  stressISPC_B(uniform int ns[],uniform  int b3,uniform int e3,
uniform int b2, uniform int e2, uniform int b1,
uniform int e1,
uniform float x[], uniform float y[],uniform float z[],
uniform float dt2,
uniform float Vx[], uniform float Vy[], uniform float Vz[],
uniform float invD[], uniform float prev[], uniform float cur[],
uniform float nex[]){

uniform int n1=ns[0],n12=ns[0]*ns[1];

uniform int i2,i3;
  for( i3=b3; i3 < e3; i3++){
    for( i2=b2; i2 < e2; i2++){
      foreach(i1=b1... e1){
        int i=i1+i2*n1+i3*n12;
         nex[i]=2*cur[i]-prev[i]+dt2*invD[i]*(
           x[4]*(Vx[i-0]-Vx[i+1])+
           x[3]*(Vx[i-1]-Vx[i+2])+
           x[2]*(Vx[i-2]-Vx[i+3])+
           x[1]*(Vx[i-3]-Vx[i+4])+
           x[0]*(Vx[i-4]-Vx[i+5])+
           y[4]*(Vy[i-0*n1]-Vy[i+0])+
           y[3]*(Vy[i-1*n1]-Vy[i+2*n1])+
           y[2]*(Vy[i-2*n1]-Vy[i+3*n1])+
           y[1]*(Vy[i-3*n1]-Vy[i+4*n1])+
           y[0]*(Vy[i-4*n1]-Vy[i+5*n1])+
           z[4]*(Vz[i-0*n12]-Vz[i+0])+
           z[3]*(Vz[i-1*n12]-Vz[i+2*n12])+
           z[2]*(Vz[i-2*n12]-Vz[i+3*n12])+
           z[1]*(Vz[i-3*n12]-Vz[i+4*n12])+
           z[0]*(Vz[i-4*n12]-Vz[i+5*n12]));
       }
    }
        
}
}

export void  velISPC_B(uniform int ns[],uniform int b3,uniform int e3,
uniform int b2, uniform int e2, uniform int b1,
uniform int e1,
uniform float x[], uniform float y[],uniform float z[],
uniform float cur[], uniform float lamda[],
uniform float Vx[], uniform float Vy[], uniform float Vz[]){

uniform int n1=ns[0],n12=ns[0]*ns[1];
uniform int i2,i3;
  for( i3=b3; i3 < e3; i3++){
    for( i2=b2; i2 < e2; i2++){
      foreach(i1=b1... e1){
        int i=i1+i2*n1+i3*n12;
        Vx[i]=lamda[i]*(
           x[4]*(cur[i-1]-cur[i+0])+
           x[3]*(cur[i-2]-cur[i+1])+
           x[2]*(cur[i-3]-cur[i+2])+
           x[1]*(cur[i-4]-cur[i+3])+
           x[0]*(cur[i-5]-cur[i+4]));

        Vy[i]=lamda[i]*(
           y[4]*(cur[i-1*n1]-cur[i+0])+
           y[3]*(cur[i-2*n1]-cur[i+1*n1])+
           y[2]*(cur[i-3*n1]-cur[i+2*n1])+
           y[1]*(cur[i-4*n1]-cur[i+3*n1])+
           y[0]*(cur[i-5*n1]-cur[i+4*n1]));

        Vz[i]=lamda[i]*(
           z[4]*(cur[i-1*n12]-cur[i+0])+
           z[3]*(cur[i-2*n12]-cur[i+1*n12])+
           z[2]*(cur[i-3*n12]-cur[i+2*n12])+
           z[1]*(cur[i-4*n12]-cur[i+3*n12])+
           z[0]*(cur[i-5*n12]-cur[i+4*n12]));
        
        }
     }
   }
}

