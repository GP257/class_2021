extern "C" {

void  stressISPC(const int nb, const int *ns, const int *b3,const int *e3, 
const int *b2, const int *e2, const int *b1, 
const int *e1,
const float *x, const float *y,const float *z, 
const float dt,
const float *Vx, const float *Vy, const float *Vz,
const float *invD, const float *prev, const float *cur,
float *next);

void  velISPC(const int nb, const int *ns, const int *b3,const int *e3, 
const int *b2, const int *e2, const int *b1, 
const int *e1,
const float *x, const float *y,const float *z, 
const float *cur, const float *lamda, 
float *Vx, float *Vy, float *Vz);



void  stressISPC_B(const int *ns, const int b3,const int e3, 
const int b2, const int e2, const int b1, 
const int e1,
const float *x, const float *y,const float *z, 
const float dt,
const float *Vx, const float *Vy, const float *Vz,
const float *invD, const float *prev, const float *cur,
float *next);

void  velISPC_B( const int *ns, const int b3,const int e3, 
const int b2, const int e2, const int b1, 
const int e1,
const float *x, const float *y,const float *z, 
const float *cur, const float *lamda, 
float *Vx, float *Vy, float *Vz);



}
