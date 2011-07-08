#include<complex.h>
#define ntl 61
#define nq 4620
#define idim 21
#define nvec 48
#define mplot 10
    
struct {
    int isy;
}system_id;

struct {
    float anu, h2m2;
    complex eye;
    float pi;
} param;

struct {
    float hbar,aindb,eryd,amee,epsb,coul,akp,ra,vcd3d,rhb,rgam,rgam1,rgam2;
    integer isck;
} par_elec;

struct {
    float hbarc, p0, pf;
} par_nucl;

struct {
    float dft, dt, dpz, dp3;
} grid_size;

struct {
       int ix[nq],jx[nq],kx[nq];
       int irq[2*idim+1,2*idim+1,2*idim+1];   
}grid_label;
       
/*struct {
    complex gre[17191020];       ( was [4620][61][61] )  gre[nq][nt1][nt1]
} green;
*/

/*struct {
    complex ftcoef[189];         ftcoef[nvec+15+nvec+15+nvec+15]
} fft_coef;
*/

/*struct {
    complex sgla[281820]        /* was [4620][61] , sgle[281820]      sgla(nq,ntl),sgle(nq,ntl)   
            was [4620][61] ;
} sigma;
*/

struct {
    int ifk;
//    complex smf[nq];
} sig_mf;

//struct {
//    float pot[79507]     /* was [43][43][43] */;       pot(-idim:idim,-idim:idim,-idim:idim)
//} pot_sqd;
/*struct {
    complex   pot_fk(nvec,nvec,nvec);
}pot_mf;*/

//struct {
//    complex dpot[nq];
//} pot_a;

struct {
       float oke(ntl),ope_corr(ntl),ope_mf(ntl),ote(ntl),oden(ntl),oqua(ntl),omega_pl;
}output;       

struct {
    int nplot, iplot[mplot];
} plot;
