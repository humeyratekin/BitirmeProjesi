#include <stdio.h>

/* var olan parametre deðerlerini ve tahmin edilen uzunluklarýiçeren tablo giriþlerini tutan yapýdýr*/

struct table_entry_structure{
	
	double u,length;
    

}table_entry_td;

/*eðrinin aralýklarýný tutan yapýdýr */
//(adaptif integrasyon fonksiyonunda kullanýlacaktýr)

typedef struct interval_structure{
	
	double u1,u2;
	double length;
}interval_td;

//2 boyutlu kübik eðriler için katsayýlarý tutar.
typedef struct cubic_curve_structure{
	double ax,bx,cx,dx;
	double ay,by,cy,dy;
}cubic_curve_td;

//polinom fonksiyon yapýsý
typedef struct polynomial_structure{
	double *coeff;
	int degree;
}polynomial_td;

/* Adaptive Integration
 eðrinin uzunluðunu hesaplayan fonksiyon
*/

void adaptive_integration(cubic_curve_td *curve,double u1, double u2,double tolerance)
{
	double subdivide();
	polynomial_td func;
	interval_td full_interval;
	double total_length;
	double integrate_func();
	double temp;
	func.degree=4;
	func.coeff=(double *)malloc(sizeof(double)*5);
	func.coeff[4]=9*(curve->ax*curve->ax+curve->ay*curve->ay);
	func.coeff[3]=12*(curve->ax*curve->bx+curve->ay*curve->by);
	func.coeff[2]=(6*(curve->ax*curve->cx+curve->ay*curve->cy)+4*(curve->bx*curve->bx+curve->by*curve->by));
	func.coeff[1]=4*(curve->bx*curve->cx+curve->by*curve->cy);
	func.coeff[0]=curve->cx*curve->cx+curve->cy*curve->cy;
	
	full_interval.u1=u1; 
	full_interval.u2=u2;
	
	temp=integrate_func(&func,&full_interval);
	printf("\nBaslangic tahmini=%1f u1: %1f u2: %1f\n",temp,u1,u2);
	full_interval.length=temp;
	total_length=subdivide(&full_interval,&func,0.0,tolerance);
	printf("\n Toplam uzunluk= %1f\n",total_length);
	
}

/*
SUBDIVIDE


*/
 

double subdivide(interval_td *full_interval,polynomial_td *func,double total_length, double tolerance)
{
	interval_td left_interval, right_interval;
	double left_length,right_length;
	double midu;
	double subdivide();
	double integrate_func();
	double temp;
	void add_table_entry();
	
	midu=(full_interval->u1+full_interval->u2)/2;
	left_interval.u1=full_interval->u1;
	left_interval.u2=midu;
	right_interval.u1=midu; 
	right_interval.u2=full_interval->u2;
	
	left_length=integrate_func(func,&left_interval);
	right_length=integrate_func(func,&right_interval);
	
	temp=fabs(full_interval->length-(left_length+right_length));
	
	if(temp>tolerance){
		left_interval.length=left_length;
		right_interval.length=right_length;
		total_length=subdivide(&left_interval,func,total_length,tolerance/2.0);
		total_length=subdivide(&right_interval,func,total_length,tolerance/2.0);
		return(total_length);

	}
	else{
		 total_length = total_length + left_length;
		 add_table_entry(midu,total_length);
		 total_length = total_length + right_length;
		 add_table_entry(full_interval->u2,total_length);
		 return(total_length);

	}
}

/* ------------------------------------------------------------------------
parametrik deðer ve yay uzunluðunu listeleyen bir tablo oluþturur
*/
void add_table_entry(double u, double length)
{
 
 printf("\ntable entry: u: %lf, length: %lf",u,length);
}

/* */

double integrate_func(polynomial_td *func,interval_td *interval)
{

double x[5]={.1488743389,.4333953942,.6794095682,.8650633666,.9739065285};
double w[5]={.2966242247,.2692667193,.2190863625,.1494513491,.0666713443};

double length,midu,dx,diff;
int i;
double evaluate_polynomial();
double u1,u2;

u1=interval->u1;
u2=interval->u2;

midu=(u1/u2)/2.0;
diff=(u2-u1)/2.0;

length=0.0;
for(i=0; i<5; i++){    
	dx=diff*x[i];
	length+=w[i]*(sqrt(evaluate_polynomial(func,midu+dx))+sqrt(evaluate_polynomial(func,midu-dx)));
}
length *=diff; //interval aralýðýnda ölçekler
return(length);
}


/*Polinom oluþturma fonksiyonu*/
double evaluate_polynomial(polynomial_td *poly, double u){
	double w;
	int i;
	double value;
	value=0.0;
	w=1.0;
	for(i=0; i<=poly->degree; i++){
		value+=poly->coeff[i]*w;
		w*=u;
	}
	return value;
	
}

int main(int argc, char** argv) {
	
	  /* Eðrinin katsayýlarý*/
    double ax, ay, bx, by, cx, cy, dx, dy;
    double u1, u2;      //Ölçülecek eðrinin intervali
    
    double tolerance=0.3;
    
    polynomial_td *x = (polynomial_td *) malloc (sizeof(polynomial_td));
    polynomial_td *y = (polynomial_td *) malloc (sizeof(polynomial_td));
    
      /*** Kübik eðriyi iki fonksiyon olarak tanýmlar ***/
    ax = 10;       //u^3 x(u) fonksiyonunun katsayýlarý
    bx = 5;         //u^2
    cx = 1;         //u
    dx = 10;       //sabit

    ay = 4;         //u^3 y(u) fonksiyonunun katsayýlarý
    by = 0;         
    cy = 1;
    dy = 0;
    
     u1 = 0.8;
    u2 = 1;
    
    cubic_curve_td *p=(cubic_curve_td*)malloc(sizeof(cubic_curve_td));
    p->ax=ax; p->ay=ay;
    p->bx=bx; p->by=by;
    p->cx=cx; p->cy=cy;
    p->dx=dx; p->dy=dy;
    
    
    
    printf("Kubik egri denklemleri:\n");
    printf("x(u) = %g * u^3 + %g * u^2 + %g * u + %g\n", ax, bx, cx, dx);
    printf("y(u) = %g * u^3 + %g * u^2 + %g * u + %g\n\n", ay, by, cy, dy);
    printf("Adaptif integrasyon:\n");
    adaptive_integration(p,u1,u2,tolerance);
    

	return 0;
}
