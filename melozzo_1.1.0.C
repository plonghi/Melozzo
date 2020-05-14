#include <passe_par_tout.h>
#define verde 2
#define bianco 7
#define cyan 5
#define rosso 1
#define giallo 4
#define blu 3
#define magenta 6
#define nero 0

using namespace std;
int menu();
int julia_mandelbrot();
int menu_mappe();
int menu_automi();
int credits();

int main()
{
	char *t="MENU";
	int K=1440, L=900;
	m_startg (t,&K,&L);
	menu();
	int B=0;
	while (B==0)
	{
	int finestra=0,m_select(finestra);
	//preparazione del disegno da effettuare
	double coord[2];
	int Xoffset=17,Yoffset[4]={160,218,275,332},Lx=121,Ly=44,Xquit=347,Yquit=20;
		m_mouse(coord);
		if (coord[0]>Xoffset&&coord[0]<Xoffset+Lx)
		{
			if (coord[1]>Yoffset[0]&&coord[1]<Yoffset[0]+Ly) cerr<<"hai scelto: credits\n",m_close(&(finestra=0)),credits(),menu();
			else if (coord[1]>Yoffset[1]&&coord[1]<Yoffset[1]+Ly)  m_close(&(finestra=0)),menu_automi(),menu();
			else if (coord[1]>Yoffset[2]&&coord[1]<Yoffset[2]+Ly) m_close(&(finestra=0)),julia_mandelbrot(),menu();
			else if (coord[1]>Yoffset[3]&&coord[1]<Yoffset[3]+Ly) m_close(&(finestra=0)),menu_mappe(),menu();
		}
		else if(coord[0]>Xquit&&coord[0]<Xquit+Lx&&coord[1]>Yquit&&coord[1]<Yquit+Ly) cerr<<"hai scelto quit!\n",B=1;
	}
	m_endg();
	return 1; 
}

int credits()
{
	char *t="credits";
	int dove[2]={400,200};
	m_place_window(dove, t);
	int K=(100 << 24) + (100 << 16) + (100 << 8) + 0,L=(100<<24)+(728<<11)+700;
	m_window (&K,&L);
	K=1,L=verde;
	m_use_as_pixmap(&K,&L);
	int opz[]={0,0,0,0,728,700,1};
	m_load_external_pixmap(&(K=1),"credits.ppm",opz);
	m_wait_for_events(&(K=999));
	m_close(&(K=0));
	return 1; 
}

int menu()
{
	//disegno il MENU
	double xmin =0,ymin=0,xmax=492,ymax=480, dove[2]={170.0,10.0};
	int K=(100 << 24) + (100 << 16) + (100 << 8) + 0,L=(100<<24)+(492<<11)+480;
	m_window (&K,&L);
	K=1,L=verde;
	m_use_as_pixmap(&K,&L);
	int opz[]={0,0,0,0,480,492,1};
	m_load_external_pixmap(&(K=1),"menu.ppm",opz);
	m_frame(&xmin, &ymin, &xmax, &ymax);
	m_color(&(K=bianco)),m_text("                             pietro.longhi@gmail.com",dove);
	//HO disegnato il menu
}

int Mandelbrot(double X,double Y,double scale,int lato, int *precisione);
int Julia(float a, float b,double X,double Y,double scale,int lato,int *precisione);

int julia_mandelbrot()
{
	//disegno il MENU
	char *t="MENU";
	const int n_pulsanti=4;
	int lato_pulsanti=100, Lx=800,Ly=800, lato_mandelbrot=300, lato_julia=300, K[3]={(100 << 24) + (100 << 16) + (100 << 8) + 100,(100 << 24) + (0 << 16) + (0 << 8) + 100, (100 << 24) + (0 << 16) + (0 << 8) + 100}, L[3]={(100<<24)+(lato_pulsanti<<11)+lato_pulsanti*n_pulsanti,(100 << 24) + (lato_mandelbrot << 11) + lato_mandelbrot, (100 << 24) + (lato_julia << 11) + lato_julia};
	double xmin =0,ymin=0,xmax=1,ymax=n_pulsanti;
	m_window (&K[0],&L[0]);
	m_frame(&xmin, &ymin, &xmax, &ymax);
	int finestra=1,colore=giallo;
	m_use_as_pixmap(&finestra,&colore);
	double vertici[8]={0,0,1,0,1,1,0,1}, coord[2]={0.4,0.4};
	int n_vertici=4;
	int puls_color[n_pulsanti]={rosso,giallo,blu,verde};
	char *menu[n_pulsanti]={"quit","choose","reset","zoom"};
	colore=nero;
	for (int p=0;p<n_pulsanti;p++)
	{
	m_color(&puls_color[p]),m_polygon(vertici, &n_vertici);
	vertici[1]++,vertici[3]++,vertici[5]++,vertici[7]++;
	m_color(&colore),m_text(menu[p],coord),coord[1]++,m_flush();
	}
	//HO disegnato il menu
	//disegno le finestre di lavoro
	t="Mandelbrot";
	int dove[2]={200,0};
	m_place_window(dove, t);
	m_window (&K[1],&L[1]);
	finestra=2,colore=nero,m_use_as_pixmap(&finestra,&colore);
	float semilato=lato_mandelbrot/2;
	xmin=-semilato,ymin=-semilato,xmax=semilato,ymax=semilato;
	m_frame(&xmin, &ymin, &xmax, &ymax);
	t="Julia";
	dove[0]=600,dove[1]=0;
	m_place_window(dove, t);
	m_window (&K[2],&L[2]);
	finestra=3,colore=nero,m_use_as_pixmap(&finestra,&colore);
	semilato=lato_julia/2;
	xmin=-semilato,ymin=-semilato,xmax=semilato,ymax=semilato;
	m_frame(&xmin, &ymin, &xmax, &ymax);
	//disegnate finestre di lavoro
	//preparazione del disegno da effettuare
	double coord_m[2]={0,0},Mscale=1.4,Jscale=2,XM=-0.8,YM=0,XJ=0,YJ=0,aJ=0,bJ=0;
	int Mprecisione=200,Jprecisione=50;
	Mandelbrot(XM,YM,Mscale,lato_mandelbrot, &Mprecisione);
	Julia(aJ,bJ,XJ,YJ,Jscale,lato_julia, &Jprecisione);
	int A=0;
	while(A==0)
	{
		m_mouse(coord_m),m_which_window(&finestra);
		if (finestra==0)
		{
			if((coord_m[1]>3)&&(coord_m[1]<4))
			{
				int B=0;
				while(B==0)
				{m_mouse(coord_m),m_which_window(&finestra);
				if(finestra==1) m_select(&finestra),XM=coord_m[0]*Mscale/(lato_mandelbrot/2)+XM,YM=coord_m[1]*Mscale/(lato_mandelbrot/2)+YM,Mscale=Mscale/2,Mandelbrot(XM,YM,Mscale,lato_mandelbrot,&Mprecisione),B=1;
				if(finestra==2) m_select(&finestra),XJ=coord_m[0]*Jscale/(lato_julia/2)+XJ,YJ=coord_m[1]*Jscale/(lato_julia/2)+YJ,Jscale=Jscale/2,Julia(aJ,bJ,XJ,YJ,Jscale,lato_julia,&Jprecisione),B=1;
				}
			};
			if((coord_m[1]>2)&&(coord_m[1]<3)) XM=-0.8,YM=XJ=YJ=aJ=bJ=0,Mscale=1.4,Jscale=2,Mprecisione=200,Jprecisione=50,Mandelbrot(XM,YM,Mscale,lato_mandelbrot, &Mprecisione), Julia(aJ,bJ,XJ,YJ,Jscale,lato_julia, &Jprecisione);
			if ((coord_m[1]>1)&&(coord_m[1]<2)) 
			{
				int B=0;
				while(B==0)
				{
					m_mouse(coord_m),m_which_window(&finestra);
					if(finestra==1) aJ=coord_m[0]*Mscale/(lato_mandelbrot/2)+XM,bJ=coord_m[1]*Mscale/(lato_mandelbrot/2)+YM,XJ=YJ=0,Jscale=2,Jprecisione=50,Julia(aJ,bJ,XJ,YJ,Jscale,lato_julia,&Jprecisione),B=1;
				}
			}					
			if ((coord_m[1]<1)&&(coord_m[1]>0)) A=1;
		}
	}
	m_close(&(finestra=2)),m_close(&(finestra=1)),m_close(&(finestra=0));
	return 1; 
}

double xprim_M(double x,double p,double a)
{
	return pow(x,2)-pow(p,2)+a;
}
double pprim_M(double x,double p,double b)
{
	return 2*x*p+b;
}

int Mandelbrot(double X,double Y,double scale,int lato, int *precisione)
{
	int max_iter=0,colore,finestra=1,s=0;       /*massime iterazioni di ciascun punto*/
	float R=1.4*sqrt(2.0);    /*ZOOM, raggio di tolleranza*/
	double x1, p1, grandezza=1,coord[2];
	m_select(&finestra);
	cerr<<"MANDELBROT---->\nPrecisione: "<<*precisione<<endl;
	cerr<<"Nuovo centro finestra:\nX= "<<X<<"\nY= "<<Y<<"\nSCALA: "<<scale<<"\nlato: "<<lato<<endl;
	m_symbol(&s),m_size_symbol(&grandezza);
		for(int i=0;i<lato;i++)
		{
			for(int j=0;j<lato;j++)
			{
				double x=0, p=0,a=(float(i)/lato)*(2*scale)-scale+X,b=(float(j)/lato)*(2*scale)-scale+Y;
				coord[0]=i-lato/2,coord[1]=j-lato/2;
				int iter=0;
				while((sqrt(pow(x,2)+pow(p,2))<R)&&(iter<*precisione))
				{
					x1=xprim_M(x,p,a);
					p1=pprim_M(x,p,b);
					x=x1;
					p=p1;
					iter++;
				};
				if ((iter>max_iter)&&(iter<*precisione)) max_iter=iter;
				if (iter==*precisione) colore=nero;
				else colore=1+int(7*(pow(float(iter)/(*precisione),0.335)));

				m_color(&colore);

				m_point(coord);
			}
		}
		if(max_iter>*precisione-10) *precisione+=100;
	return 1;
}

int Julia(float a, float b,double X,double Y,double scale,int lato,int *precisione)
{
	int max_iter=0,iter=0,colore=nero,finestra=2,s=0;       /*massime iterazioni di ciascun punto*/
	double x1, p1,coord[2],R=15,grandezza=1;
	m_select(&finestra);
	cerr<<"JULIA---->\nPrecisione: "<<*precisione<<"\nparametro: "<<a<<" + i"<<b<<endl;
	cerr<<"Nuovo centro finestra:\nX= "<<X<<"\nY= "<<Y<<"\nSCALA: "<<scale<<"\nlato: "<<lato<<endl;
	m_symbol(&s),m_size_symbol(&grandezza);
	for(int i=0;i<lato;i++)
	{
		for(int j=0;j<lato;j++)
		{
			double x=(float(i)/lato)*(2*scale)-scale+X, p=(float(j)/lato)*(2*scale)-scale+Y;
			coord[0]=i-lato/2,coord[1]=j-lato/2;
			iter=0;
			for(int t=0;t<*precisione;t++)
			{
				x1=xprim_M(x,p,a);
				p1=pprim_M(x,p,b);
				x=x1;
				p=p1;
				iter++;
				if (sqrt(pow(x,2)+pow(p,2))>R) break;
			};
			if ((iter>max_iter)&&(iter<*precisione)) max_iter=iter;
			if (iter==*precisione) colore=nero;
			else colore=1+int(7*(pow(float(iter)/(*precisione),0.55)));
			m_color(&colore);
			m_point(coord);
		}
	}
	if(max_iter>*precisione-20) *precisione+=100;
	return 0;
}

int mappe();
int mappa_generalizzata();
int mappe_a(int n,int iter,double scale,double betha, int finestra,double X,double Y, int K,int L,int k);
int mappe_b(int n,int iter,double scale,double alpha,double betha,int finestra,double X,double Y, int K,int L);
int aggiungi_punti(int N,double scale,double X,double Y);
double xprim_a(double x,double y,double alpha);
double pprim_a(double x,double y,double alpha);
double xprim_b(double x,double y,double a);
double pprim_b(double x,double y,double b);
double xinv_a(double x,double y,double alpha);
double pinv_a(double x,double y,double alpha);
double perturb(double x,double y);

int menu_mappe()
{
	//disegno il MENU
	char *t="MENU";
	const int n_pulsanti=3;
	int lato_pulsanti=100, K={(100 << 24) + (100 << 16) + (100 << 8) + 100}, L={(100<<24)+(2*lato_pulsanti<<11)+lato_pulsanti*n_pulsanti};
	double xmin =0,ymin=0,xmax=1,ymax=n_pulsanti;
	m_window (&K,&L);
	m_frame(&xmin, &ymin, &xmax, &ymax);
	int finestra=1,colore=giallo;
	m_use_as_pixmap(&finestra,&colore);
	double vertici[8]={0,0,1,0,1,1,0,1}, coord[2]={0.2,0.4};
	int n_vertici=4;
	int puls_color[n_pulsanti]={giallo,blu,verde};
	char *menu[n_pulsanti]={"menu principale","mappa quadratica","generalizzata quadratica"};
	colore=nero;
	for (int p=0;p<n_pulsanti;p++)
	{
		m_color(&puls_color[p]),m_polygon(vertici, &n_vertici);
		vertici[1]++,vertici[3]++,vertici[5]++,vertici[7]++;
		m_color(&colore),m_text(menu[p],coord),coord[1]++,m_flush();
	}
	//disegnato MENU
	double coord_m[2]={0,0};
	int A=0;
	while(A==0)
	{
		m_mouse(coord_m),m_which_window(&finestra);
		if (finestra==0)
		{
			if ((coord_m[1]>2)&&(coord_m[1]<3)){m_close(&(finestra=0)),mappa_generalizzata(),A=1;} 
			else if ((coord_m[1]>1)&&(coord_m[1]<2)){m_close(&(finestra=0)),mappe(),A=1;}
			else if ((coord_m[1]<1)&&(coord_m[1]>0)) {m_close(&(finestra=0)),A=1;}
		}
	}
	return 1; 
}
int mappe()
{
	//disegno il MENU
	double scale=3,betha=0.205,X=0,Y=0;
	int iter,N_punti=0;
	char *t="MENU";
	const int n_pulsanti=10;
	int lato_x_pulsanti=100, lato_y_pulsanti=50, Lx=600,Ly=600, K[3]={(100 << 24) + (100 << 16) + (100 << 8) + 100,(100 << 24) + (0 << 16) + (0 << 8) + 100,(100 << 24) + (0 << 16) + (180 << 8)}, L[3]={(100<<24)+(lato_x_pulsanti<<11)+lato_y_pulsanti*n_pulsanti,(100 << 24) + (Lx << 11) + Ly,(100<<24)+(400<<11)+70}, s=10;
	double xmin =0,ymin=0,xmax=1,ymax=n_pulsanti;
	m_window (&K[0],&L[0]);
	m_frame(&xmin, &ymin, &xmax, &ymax);
	int finestra=1,colore=giallo;
	m_use_as_pixmap(&finestra,&colore);
	double vertici[8]={0,0,1,0,1,1,0,1}, coord[2]={0.15,0.4},coord_m[2];
	int n_vertici=4;
	int puls_color[n_pulsanti]={rosso,giallo,blu,verde,magenta,cyan,giallo,blu,verde,magenta};
	char *menu[n_pulsanti]={"quit","reset","clear","zoom out","chs.pnt","zoom","autofill","iter","alfa","draw G(k)"};
	colore=nero;
	for (int p=0;p<n_pulsanti;p++)
	{
		m_color(&puls_color[p]),m_polygon(vertici, &n_vertici);
		vertici[1]++,vertici[3]++,vertici[5]++,vertici[7]++;
		m_color(&colore),m_text(menu[p],coord),coord[1]++,m_flush();
	};
	//HO disegnato il menu
	//disegno menu alfa
	cerr<<"\nalfa vale di default 0.205\n",t="Alfa";
	int dove[2]={1000,200};
	xmin =-0.2,ymin=-0.5,xmax=0.7,ymax=0.5;
	m_place_window(dove, t);
	m_window (&K[2],&L[2]);
	m_frame(&xmin, &ymin, &xmax, &ymax), m_color(&(L[3]=bianco)), t="0.5", coord[0]=0.515,coord[1]=0, m_text(t, coord), t="+",coord[0]=0.6, m_text(t, coord), t="-",coord[0]=-0.1, m_text(t, coord), t="0", coord[0]=-0.03, m_text(t, coord), coord[0]=coord[1]=0, m_move(coord), coord[0]=0.5, m_thickness(&(K[3]=4)),m_line(coord);
	//HO disegnato il menu alfa
	//disegno finestra della mappa
	t="Mappa";
	dove[0]=200,dove[1]=0;
	m_place_window(dove, t);
	m_window (&K[1],&L[1]);
	colore=bianco,finestra=2,finestra=3,m_use_as_pixmap(&finestra,&colore);
	xmin=-scale,ymin=-scale,xmax=scale,ymax=scale;
	m_frame(&xmin, &ymin, &xmax, &ymax);
	//HO disegnato finestra della mappa
	FILE* file=fopen ("punti_partenza","wt");
	fclose(file);
	int A=0,G=1;
	iter=10000;
	while(A==0)
	{
		system("clear");
		cout<<"Orbite: "<<N_punti<<"\nIterazioni: "<<iter<<"\nScala: "<<scale<<"\nCentro finestra: ("<<X<<","<<Y<<")\nValore di alfa: "<<betha<<"\nDisegna ogni: "<<G<<" iterazioni\n";
		m_mouse(coord_m),m_which_window(&finestra);
		if (finestra==0)
		{
			if ((coord_m[1]<10)&&(coord_m[1]>9))
			{
				cerr<<"Quale iterata vuoi vedere? ";
				int control=0;
				while(control==0)
				{
					cin>>G;
					if(betha>0) control=1;
					else cerr<<"\nErrore, devi partire almeno dalla 1a iterata!! ";
				}
				finestra=2,m_select(&finestra),m_close(&finestra),mappe_a(N_punti,iter,scale,betha,finestra,X,Y,K[1],L[1],G);
			};//scegli quale iterata vuoi vedere
			if ((coord_m[1]<9)&&(coord_m[1]>8))
			{
				cerr<<"hai scelto immissione di alfa\nImmettere alfa (tra 0 e 1/2): ";
				int control=0;
				while(control==0)
				{
				cin>>betha;
				if((betha>=0)&&(betha <=0.5)) control=1;
				else cerr<<"\nErrore, inserire un opportuno valore di alfa compreso tra 0 e 1/2!! ";
				}
				finestra=2,m_select(&finestra),m_close(&finestra),mappe_a(N_punti,iter,scale,betha,finestra,X,Y,K[1],L[1],G);
			};//immetti alfa
			if ((coord_m[1]<8)&&(coord_m[1]>7))
			{
				cerr<<"hai scelto immissione di iterazioni\nImmettere quante iterazioni: ";
				int control=0;
				while(control==0)
				{
				cin>>iter;
				if(iter>0)control=1;
				else cerr<<"\nErrore, inserire un opportuno valore maggiore di zero per le iterate!!! ";
				}
				finestra=2,m_select(&finestra),m_close(&finestra),mappe_a(N_punti,iter,scale,betha,finestra,X,Y,K[1],L[1],G);
			};//iter
			if ((coord_m[1]<7)&&(coord_m[1]>6))
			{
				cerr<<"hai scelto autofill\n";
				finestra=2,aggiungi_punti(300,scale,X,Y),m_select(&finestra),m_close(&finestra),mappe_a((N_punti+=300),iter,scale,betha,finestra,X,Y,K[1],L[1],G);
			};//autofill
			if ((coord_m[1]<6)&&(coord_m[1]>5))
			{
				cerr<<"hai scelto zoom\n";
				int B=0;
				while(B==0)
				{m_mouse(coord),m_which_window(&finestra);
				if(finestra==2) m_select(&finestra),X=coord[0],Y=coord[1],scale=scale/2,m_close(&finestra),mappe_a(N_punti,iter,scale,betha,finestra,X,Y,K[1],L[1],G);
				if(finestra==0||finestra==1) B=1;
				};
			};//zoom
			if ((coord_m[1]<5)&&(coord_m[1]>4))
			{
				cerr<<"hai scelto choose point\n";
				int B=0;
				while(B==0)
				{m_mouse(coord),m_which_window(&finestra);
				if(finestra==2)
				{
					double a,b,s=1,c=coord[0],d=coord[1];
					FILE* file=fopen ("punti_partenza","a");
					fprintf(file,"%e        %e\n",c,d),N_punti+=1;
					fclose(file);
					m_symbol(&(B=0)),m_size_symbol(&s),m_color(&(B=nero)),m_select(&finestra),m_point(coord);
					for(int j=0;j<iter/2;j++)
					{
						for(int y=0;y<G;y++)
						{
						a=xprim_a(coord[0],coord[1],2*M_PI*betha);
						b=pprim_a(coord[0],coord[1],2*M_PI*betha);
						coord[0]=a;
						coord[1]=b;
						}
						if(sqrt(pow(coord[0]-X,2)+pow(coord[1]-Y,2))<100) m_point(coord);
						else break;
					};
					coord[0]=c,coord[1]=d;
					for(int j=0;j<iter/2;j++)
					{
						for(int y=0;y<G;y++)
						{
						a=xinv_a(coord[0],coord[1],2*M_PI*betha);
						b=pinv_a(coord[0],coord[1],2*M_PI*betha);
						coord[0]=a;
						coord[1]=b;
						}
						if(sqrt(pow(coord[0]-X,2)+pow(coord[1]-Y,2))<100) m_point(coord);
						else break;
					};
					B=0;
				}
				if (finestra!=2) B=1;
				};
			};//choose pnt
			if((coord_m[1]>3)&&(coord_m[1]<4))
			{
				cerr<<"hai scelto zoom out\n";
				finestra=2,m_select(&finestra),scale=scale*2,m_close(&finestra),mappe_a(N_punti,iter,scale,betha,finestra,X,Y,K[1],L[1],G);
			};//zoom out
			if((coord_m[1]>2)&&(coord_m[1]<3))
			{
				FILE* file=fopen ("punti_partenza","wt");
				fclose(file);
				N_punti=0,finestra=2,m_close(&finestra),mappe_a(N_punti,iter,scale,betha,finestra,X,Y,K[1],L[1],G),finestra=0;
			};//clear
			if ((coord_m[1]>1)&&(coord_m[1]<2))
			{
				FILE* file=fopen ("punti_partenza","wt");
				fclose(file);
				N_punti=0,X=Y=0,scale=3,iter=10000,G=1,finestra=2,m_close(&finestra),mappe_a(N_punti,iter,scale,betha,finestra,X,Y,K[1],L[1],G),finestra=0;
			}; //reset	
			if ((coord_m[1]<1)&&(coord_m[1]>0)) A=1;
		}
		else if (finestra==1)
		{
			if ((coord_m[0]<=0.5)&&(coord_m[0]>=0)&&(coord_m[1]<=0.2)&&(coord_m[1]>=-0.2)) betha=coord_m[0],cerr<<"alfa vale "<<betha<<endl;
			else if ((coord_m[0]<=-0.05)&&(coord_m[0]>=-0.15)&&(coord_m[1]<=0.2)&&(coord_m[1]>=-0.2)) 
			{
				if(betha>0.001) betha-=0.00005;
				cerr<<"decremento:\talfa vale "<<betha<<endl;
			}
			else if ((coord_m[0]<=0.65)&&(coord_m[0]>=0.55)&&(coord_m[1]<=0.2)&&(coord_m[1]>=-0.2)) 
			{
				if(betha<0.499) betha+=0.00005;
				cerr<<"aumento:\talfa vale "<<betha<<endl;
			}
			m_close(&(finestra=2)),mappe_a(N_punti,iter,scale,betha,finestra,X,Y,K[1],L[1],G),finestra=0;
		};
	};
	m_close(&(finestra=0)),m_close(&(finestra=1)),m_close(&(finestra=2));
	//m_endg();
	return 1; 
}


int aggiungi_punti(int N,double scale,double X,double Y) //questa funzione è fatta per aggiungere punti al file punti_partenza delle mappe, per vedere i punti, ricordarsi di chiamare la funzione mappe con argomento non più "n"ma "n+N"!!!
{
	FILE* file=fopen ("punti_partenza","a");
	for (int i=0;i<N;i++)
	{
		double x=(double) rand();
		double p=(double) rand();
		x = 2*scale*(x / 0x7fffffff) - scale + X;
		p = 2*scale*(p / 0x7fffffff) - scale + Y;
		fprintf(file,"%e        %e\n",x,p);
	};
	fclose(file);
}

double perturb(double x,double y)
{
	return x*x;
}
double xprim_a(double x,double y,double alpha)
{
	return x*cos(alpha)+(y+perturb(x,y))*sin(alpha);
}
double pprim_a(double x,double y,double alpha)
{
	return (-x*sin(alpha)+(y+perturb(x,y))*cos(alpha));
}
double xinv_a(double x,double y,double alpha)
{
	return x*cos(alpha)-y*sin(alpha);
}
double pinv_a(double x,double y,double alpha)
{
	return x*sin(alpha)+y*cos(alpha)-perturb(xinv_a(x,y,alpha),y);
}

int mappe_a(int n,int iter,double scale,double betha, int finestra,double X,double Y, int K,int L,int k)
{
	char *t="Mappa";
	int dove[2]={200,0};
	m_place_window(dove, t);
	m_window (&K,&L);
	int colore=bianco;
	finestra+=1,
	m_use_as_pixmap(&finestra,&colore);
	double xmin=-scale+X,ymin=-scale+Y,xmax=scale+X,ymax=scale+Y;
	m_frame(&xmin, &ymin, &xmax, &ymax);
	double s=1;
	m_symbol(&(K=0)),m_size_symbol(&s),m_color(&(K=nero));
	double alpha=betha*2*3.141592654;
	ifstream leggi ("punti_partenza");
	for (int i=0;i<n;i++)
	{
		double a, b, c, d, coord[2];
		leggi>>coord[0];
		leggi>>coord[1];
		c=coord[0],d=coord[1];	
		m_point(coord);
		for(int j=0;j<iter/2;j++)
		{
			for(int y=0;y<k;y++)
			{
			a=xprim_a(coord[0],coord[1],alpha);
			b=pprim_a(coord[0],coord[1],alpha);
			coord[0]=a;
			coord[1]=b;
			}
			if(sqrt(pow(coord[0]-X,2)+pow(coord[1]-Y,2))<100) m_point(coord);
			else break;
		};
		coord[0]=c,coord[1]=d;
		for(int j=0;j<iter/2;j++)
		{
			for(int y=0;y<k;y++)
			{
			a=xinv_a(coord[0],coord[1],alpha);
			b=pinv_a(coord[0],coord[1],alpha);
			coord[0]=a;
			coord[1]=b;
			}
			if (sqrt(pow(coord[0]-X,2)+pow(coord[1]-Y,2))<100) m_point(coord);
			else break;
		};
	};
	leggi.close();
	return 0;
}

int mappa_generalizzata()
{
	//disegno il MENU
	double scale=3,betha=0.205,X=0,Y=0;
	int iter=10000,N_punti=0;
	char *t="MENU";
	const int n_pulsanti=9;
	int lato_x_pulsanti=100, lato_y_pulsanti=50, Lx=600,Ly=600, K[3]={(100 << 24) + (100 << 16) + (100 << 8) + 100,(100 << 24) + (0 << 16) + (0 << 8) + 100,(100 << 24) + (0 << 16) + (180 << 8)}, L[3]={(100<<24)+(lato_x_pulsanti<<11)+lato_y_pulsanti*n_pulsanti,(100 << 24) + (Lx << 11) + Ly,(100<<24)+(400<<11)+70}, s=10;
	double xmin =0,ymin=0,xmax=1,ymax=n_pulsanti;
	m_window (&K[0],&L[0]);
	m_frame(&xmin, &ymin, &xmax, &ymax);
	int finestra=1,colore=giallo;
	m_use_as_pixmap(&finestra,&colore);
	double vertici[8]={0,0,1,0,1,1,0,1}, coord[2]={0.15,0.4},coord_m[2];
	int n_vertici=4;
	int puls_color[n_pulsanti]={rosso,giallo,blu,verde,magenta,cyan,giallo,blu,verde};
	char *menu[n_pulsanti]={"quit","reset","clear","zoom out","chs.pnt","zoom","autofill","iter","a,b"};
	colore=nero;
	for (int p=0;p<n_pulsanti;p++)
	{
		m_color(&puls_color[p]),m_polygon(vertici, &n_vertici);
		vertici[1]++,vertici[3]++,vertici[5]++,vertici[7]++;
		m_color(&colore),m_text(menu[p],coord),coord[1]++,m_flush();
	};
	//HO disegnato il menu
	//disegno finestra della mappa
	t="Mappa";
	int dove[2]={200,0};
	m_place_window(dove, t);
	m_window (&K[1],&L[1]);
	colore=bianco,finestra=2,m_use_as_pixmap(&finestra,&colore);
	xmin=-scale,ymin=-scale,xmax=scale,ymax=scale;
	m_frame(&xmin, &ymin, &xmax, &ymax);
	//HO disegnato finestra della mappa
	FILE* file=fopen ("punti_partenza","wt");
	fclose(file);
	int A=0;
	float alfa=-1.4,beta=0.3;
	while(A==0)
	{
		system("clear");
		cout<<"Orbite: "<<N_punti<<"\nIterazioni: "<<iter<<"\nScala: "<<scale<<"\nCentro finestra: ("<<X<<","<<Y<<")\nValore di alfa: "<<betha<<"\nAlfa: "<<alfa<<" beta: "<<beta<<endl;
		m_mouse(coord_m),m_which_window(&finestra);
		if (finestra==0)
		{
			if ((coord_m[1]<9)&&(coord_m[1]>8))
			{
				cerr<<"hai scelto immissione di parametri\nImmettere a: ",cin>>alfa;
				cerr<<"Immettere b: ",cin>>beta;
				finestra=1,m_select(&finestra),m_close(&finestra),mappe_b(N_punti,iter,scale,alfa,beta,finestra,X,Y,K[1],L[1]);
			};//immetti alfa
			if ((coord_m[1]<8)&&(coord_m[1]>7))
			{
				cerr<<"hai scelto immissione di iterazioni\nImmettere quante iterazioni: ";
				int control=0;
				while(control==0)
				{
					cin>>iter;
					if(iter>0)control=1;
					else cerr<<"\nErrore, inserire un opportuno valore maggiore di zero per le iterate!!! ";
				}
				finestra=1,m_select(&finestra),m_close(&finestra),mappe_b(N_punti,iter,scale,alfa,beta,finestra,X,Y,K[1],L[1]),cerr<<"Iterazioni: "<<iter<<endl;
			};//iter
			if ((coord_m[1]<7)&&(coord_m[1]>6))
			{
				cerr<<"hai scelto autofill\n";
				finestra=1,aggiungi_punti(300,scale,X,Y),m_select(&finestra),m_close(&finestra),mappe_b((N_punti+=300),iter,scale,alfa,beta,finestra,X,Y,K[1],L[1]);
			};//autofill
			if ((coord_m[1]<6)&&(coord_m[1]>5))
			{
				cerr<<"hai scelto zoom\n";
				int B=0;
				while(B==0)
				{m_mouse(coord),m_which_window(&finestra);
				if(finestra==1) m_select(&finestra),X=coord[0],Y=coord[1],scale=scale/2,m_close(&finestra),mappe_b(N_punti,iter,scale,alfa,beta,finestra,X,Y,K[1],L[1]);
				if(finestra==0) B=1;
				};
			};//zoom
			if ((coord_m[1]<5)&&(coord_m[1]>4))
			{
				cerr<<"hai scelto choose point\n";
				int B=0;
				while(B==0)
				{m_mouse(coord),m_which_window(&finestra);
				if(finestra==1)
				{
					double a,b,s=1,c=coord[0],d=coord[1];
					FILE* file=fopen ("punti_partenza","a");
					fprintf(file,"%e        %e\n",c,d),N_punti+=1;
					fclose(file);
					m_symbol(&(B=0)),m_size_symbol(&s),m_color(&(B=nero)),m_select(&finestra),m_point(coord);
					for(int j=0;j<iter;j++)
					{
						a=xprim_b(coord[0],coord[1],alfa);
						b=pprim_b(coord[0],coord[1],beta);
						coord[0]=a;
						coord[1]=b;
						if(sqrt(pow(coord[0]-X,2)+pow(coord[1]-Y,2))<100) m_point(coord);
						else break;
					};
					B=0;
				}
				if (finestra!=1) B=1;
				};
			};//choose pnt
			if((coord_m[1]>3)&&(coord_m[1]<4))
			{
				cerr<<"hai scelto zoom out\n";
				finestra=1,m_select(&finestra),scale=scale*2,m_close(&finestra),mappe_b(N_punti,iter,scale,alfa,beta,finestra,X,Y,K[1],L[1]);
			};//zoom out
			if((coord_m[1]>2)&&(coord_m[1]<3))
			{
				FILE* file=fopen ("punti_partenza","wt");
				fclose(file);
				N_punti=0,finestra=1,m_close(&finestra),mappe_b(N_punti,iter,scale,alfa,beta,finestra,X,Y,K[1],L[1]),finestra=0;
			};//clear
			if ((coord_m[1]>1)&&(coord_m[1]<2))
			{
				FILE* file=fopen ("punti_partenza","wt");
				fclose(file);
				N_punti=0,X=Y=0,scale=3,iter=10000,finestra=1,m_close(&finestra),mappe_b(N_punti,iter,scale,alfa,beta,finestra,X,Y,K[1],L[1]),finestra=0;
			}; //reset	
			if ((coord_m[1]<1)&&(coord_m[1]>0)) A=1;
		}
	};
	m_close(&(finestra=0)),m_close(&(finestra=1));
	return 1; 
}

double xprim_b(double x,double y,double a)
{
	return 1-y+a*x*x;
}
double pprim_b(double x,double y,double b)
{
	return b*x;
}

int mappe_b(int n,int iter,double scale,double alpha,double betha, int finestra,double X,double Y, int K,int L) 
{
	char *t="Mappa";
	int dove[2]={200,0};
	m_place_window(dove, t);
	m_window (&K,&L);
	int colore=bianco;
	finestra+=1,
	m_use_as_pixmap(&finestra,&colore);
	double xmin=-scale+X,ymin=-scale+Y,xmax=scale+X,ymax=scale+Y;
	m_frame(&xmin, &ymin, &xmax, &ymax);
	double s=1;
	m_symbol(&(K=0)),m_size_symbol(&s),m_color(&(K=nero));
	ifstream leggi ("punti_partenza");
	for (int i=0;i<n;i++)
	{
		double a, b, c, d, coord[2];
		leggi>>coord[0];
		leggi>>coord[1];
		c=coord[0],d=coord[1];	
		m_point(coord);
		for(int j=0;j<iter/2;j++)
		{
			a=xprim_b(coord[0],coord[1],alpha);
			b=pprim_b(coord[0],coord[1],betha);
			coord[0]=a;
			coord[1]=b;
			if(sqrt(pow(coord[0]-X,2)+pow(coord[1]-Y,2))<100) m_point(coord);
			else break;
		};
	};
	leggi.close();
	return 0;
}

int automa_malattia();
int automa();
int menu_automi()
{
	//disegno il MENU
	char *t="MENU";
	const int n_pulsanti=3;
	int lato_pulsanti=100, K={(100 << 24) + (100 << 16) + (100 << 8) + 100}, L={(100<<24)+(3*lato_pulsanti<<11)+lato_pulsanti*n_pulsanti};
	double xmin =0,ymin=0,xmax=1,ymax=n_pulsanti;
	m_window (&K,&L);
	m_frame(&xmin, &ymin, &xmax, &ymax);
	int finestra=1,colore=giallo;
	m_use_as_pixmap(&finestra,&colore);
	double vertici[8]={0,0,1,0,1,1,0,1}, coord[2]={0.4,0.4};
	int n_vertici=4;
	int puls_color[n_pulsanti]={giallo,blu,verde};
	char *menu[n_pulsanti]={"menu principale","automa binario","epidemia probabilistica"};
	colore=nero;
	for (int p=0;p<n_pulsanti;p++)
	{
		m_color(&puls_color[p]),m_polygon(vertici, &n_vertici);
		vertici[1]++,vertici[3]++,vertici[5]++,vertici[7]++;
		m_color(&colore),m_text(menu[p],coord),coord[1]++,m_flush();
	}
	//disegnato MENU
	double coord_m[2]={0,0};
	int A=0;
	while(A==0)
	{
		m_mouse(coord_m),m_which_window(&finestra);
		if (finestra==0)
		{
			if ((coord_m[1]>2)&&(coord_m[1]<3)){m_close(&(finestra=0)),automa_malattia(),A=1;} 
			else if ((coord_m[1]>1)&&(coord_m[1]<2)){m_close(&(finestra=0)),automa(),A=1;}
			else if ((coord_m[1]<1)&&(coord_m[1]>0)){m_close(&(finestra=0)),A=1;}
		}
	}
	return 1; 
}

unsigned int howmany(int a, int b);
int maskbuild(unsigned int *m,int n);
int automa_0(int N,int iter,int regola);

int automa()
{
	system("clear");
	cout<<"\n------------------------------\n\t\tAUTOMA CELLULARE\n------------------------------\nquesta funzione genera un automa cellulare binario toroidale, la lunghezza d'interazione è 1 e la dinamica è totalista\n------------------------------\n";
	//disegno il MENU
	double scale=3,betha=0.205,X=0,Y=0;
	char *t="MENU";
	const int n_pulsanti=5;
	int lato_x_pulsanti=100, lato_y_pulsanti=50,K=(100 << 24) + (100 << 16) + (100 << 8) + 100, L=(100<<24)+(lato_x_pulsanti<<11)+lato_y_pulsanti*n_pulsanti, s=10;
	double xmin =0,ymin=0,xmax=1,ymax=n_pulsanti;
	m_window (&K,&L);
	m_frame(&xmin, &ymin, &xmax, &ymax);
	int finestra=1,colore=giallo;
	m_use_as_pixmap(&finestra,&colore);
	double vertici[8]={0,0,1,0,1,1,0,1}, coord[2]={0.15,0.4},coord_m[2];
	int n_vertici=4;
	int puls_color[n_pulsanti]={rosso,giallo,blu,verde,magenta};
	char *menu[n_pulsanti]={"quit","autogenerate","seed","N,iter","evolution rule"};
	colore=nero;
	for (int p=0;p<n_pulsanti;p++)
	{
		m_color(&puls_color[p]),m_polygon(vertici, &n_vertici);
		vertici[1]++,vertici[3]++,vertici[5]++,vertici[7]++;
		m_color(&colore),m_text(menu[p],coord),coord[1]++,m_flush();
	};
	//HO disegnato il menu
	int N=250,iter=200,A=0,regola=6;
	double seed=0;
	srand(time(0)),automa_0(N,iter,regola);
	while(A==0)
	{
		m_mouse(coord_m),m_which_window(&finestra);
		if (finestra==0)
		{
			if((coord_m[1]>4)&&(coord_m[1]<5))
			{
				int control=0;
				while (control==0)
				{
					cout<<"Che regola vuoi usare? (intero tra 0 e 15) ",cin>>regola;
					if (regola<0||regola>15) {system("clear"),cerr<<"la regola deve essere un numero compreso tra 0 e 15!!"<<endl,control=0;}
					else control =1;
				};
				finestra=1,m_close(&finestra),automa_0(N,iter,regola);
			};//N,iter
			if((coord_m[1]>3)&&(coord_m[1]<4))
			{
				int control=0;
				while (control==0)
				{
					cout<<"Quanto vuoi che sia lungo l'automa? ",cin>>N;
					if (N<3) {system("clear"),cerr<<"N deve valere almeno 3!!"<<endl,control=0;}
					else control =1;
				};
				cout<<"Quante volte vuoi iterare l'automa? ",cin>>iter;
				finestra=1,m_close(&finestra),automa_0(N,iter,regola);
			};//N,iter
			if((coord_m[1]>2)&&(coord_m[1]<3))
			{
				cout<<"Immetti il seed: ",cin>>seed,srand(int(seed)),finestra=1,m_close(&finestra),automa_0(N,iter,regola);
			};//seed
			if ((coord_m[1]>1)&&(coord_m[1]<2))
			{
				srand(time(0)),finestra=1,m_close(&finestra),automa_0(N,iter,regola);
			}; //autogenerate
			if ((coord_m[1]<1)&&(coord_m[1]>0)) A=1;
		}
	};
	
	m_close(&(finestra=0)),m_close(&(finestra=1));
	return EXIT_SUCCESS;
}

int maskbuild(unsigned int *m,int n)
{
	m[1]=7;
	for(int i=2;i<n;i++)//ciclo costruzione ricorsiva maschere
	{
		m[i]=m[i-1]<<1;
	};
}

unsigned int howmany(int a, int b)
{
	if (a%b==0) {return ((int)(a/b));}
	else {return ((int)(a/b)+1);}
}
int plot_automa(int N, int iter)
{
	char *t="Automa cellulare";
	ifstream file ("automa");
	char automa;
	double xmin =-double(N),ymin=-double(iter+1),xmax=0.0,ymax=0.0, grandezza=1;
	int K=(100 << 24) + (100 << 16) + (100 << 8) + 100,s=1,L=(100<<24)+(3*N<<11)+3*(iter+1);
	m_window (&K,&L);
	int finestra =2;
	int colore=cyan;
	m_use_as_pixmap (&finestra,&colore);
	m_frame(&xmin, &ymin, &xmax, &ymax);
	double coord[2];
	m_size_symbol(&grandezza);
	colore=nero;
	m_color(&colore);
	m_symbol(&s);
	for(int y=0;y>-(iter+1);y--)
	{for(int x=0;x>-N;x--)
	{
		coord[0]=x,coord[1]=y;
		file>>automa; 
		m_move(coord);
		if(automa=='1'){colore=verde,m_color(&colore),m_point(coord);}
		else if(automa=='0'){colore=bianco,m_color(&colore),m_point(coord);}
		else {colore=magenta,m_color(&colore),m_point(coord);};
	}
	}
	file.close();
	//m_wait_for_events(&(L=0));
	//m_close(&(finestra=1));
	return 1; 
}

int automa_0(int N,int iter,int regola)
{
	ofstream file("automa"); 
	unsigned int k=howmany(N,32);
	unsigned int *mask=new unsigned int[N];
	maskbuild(mask,N);
	unsigned int *a=new unsigned int[k];
	for (int i=0;i<k;i++) 
	{
		if (i==(k-1)) a[i]=(int)(((double)rand()/RAND_MAX)*((int)pow(2.0,N%32)-1));
		else a[i]=(int)(((double)rand()/RAND_MAX)*((int)pow(2.0,32)-1));
	};
	unsigned int bit[k*N];
	for(int z=k-1;z>=0;z--)//stampa binaria di a
	{
		int A=32;   
		if(z==k-1) A=N%32;    
		for(int i=A-1;i>=0;i--)
		{
			unsigned int c=((unsigned int)pow(2.0,i)&a[z])>>i;
			file<<c;
		}
	}
	file<<endl<<endl;
	int R[4];
	cout<<"La regola è la "<<regola<<" ovvero ";
	for (int i=3;i>=0;i--)
	{
		R[i]=(int)(((float)regola)/pow(2.0,i)),cout<<R[i];
		if (R[i]!=0) regola-=R[i]*(int)pow(2.0,i);
	}
	cout<<endl;
	for (int j=0;j<iter;j++)
	{
		unsigned int a1[k];  
		int A[k];//lunghezze delle stringhe k-1esime
		for(int l=0;l<k-1;l++) 
		{
			A[l]=32;
		};       
		A[k-1]=N%32;    
		for(int z=k-1;z>=0;z--)
		{
			a1[z]=0;
			unsigned int c=0;
			for(int i=(A[z]-1);i>=0;i--)//confronto con maschere e scrittura del successivo automa!!
			{
				if(i>=1&&i<=A[z]-2) c=(mask[i]&a[z])>>i-1;//caso normale
				else if (i==A[z]-1) {if (z==(k-1)) c=(((unsigned int)(pow(2.0,A[k-1]-1)+pow(2.0,A[k-1]-2))&a[k-1])>>A[k-1]-2)+((1&a[0])<<2);//ultima stringa, ultimi numeri
					else c=(((unsigned int)(pow(2.0,A[z]-1)+pow(2.0,A[z]-2))&a[z])>>A[z]-2)+((1&a[z+1])<<2);//generica stringa, ultimi numeri
				}
				else if (i==0) {if (z==(0)) c=(((unsigned int)pow(2.0,A[k-1]-1)&a[k-1])>>A[k-1]-1)+((3&a[0])<<1);
					else c=(((unsigned int)pow(2.0,A[z-1]-1)&a[z-1])>>A[z-1]-1)+((3&a[z])<<1);
				}
				else cerr<<"non va bene in automa cellulare!!"<<endl,exit(255);
				switch(c){
					case 7: bit[z*32+i]=R[3];
					break;
					case 5: case 3: case 6: bit[z*32+i]=R[2];
					break;
					case 1:case 2:case 4: bit[z*32+i]=R[1];
					break;
					case 0: bit[z*32+i]=R[0];
					break;
					default: bit[z*32+i]=8;
					break;}
					file<<bit[z*32+i];
					a1[z]=a1[z]+bit[z*32+i]*(int)pow(2.0,i);
			}
		};
		for(int i=0;i<k;i++) {a[i]=a1[i];};
		file<<endl;
	}
	delete [] mask;
	delete [] a;
	file.close(),
	plot_automa(N,iter);
}
int plot_malattia (int N,int iter,int lungint,float PG,float PC,float PD,float PN);
int automa_malattia_0(int L,string a,float PG,float PC,float PD,float PN);
int automa_malattia()
{
	system("clear");
	cout<<"\n------------------------------\n\t\tAUTOMA CELLULARE\n------------------------------\nquesta funzione genera un automa cellulare probabilistico toroidale, la dinamica è totalista\n------------------------------\n";
	//disegno il MENU
	double scale=3,betha=0.205,X=0,Y=0;
	char *t="MENU";
	const int n_pulsanti=4;
	int lato_x_pulsanti=150, lato_y_pulsanti=50,K=(100 << 24) + (100 << 16) + (100 << 8) + 100, L=(100<<24)+(lato_x_pulsanti<<11)+lato_y_pulsanti*n_pulsanti, s=10;
	double xmin =0,ymin=0,xmax=1,ymax=n_pulsanti;
	m_window (&K,&L);
	m_frame(&xmin, &ymin, &xmax, &ymax);
	int finestra=1,colore=giallo;
	m_use_as_pixmap(&finestra,&colore);
	double vertici[8]={0,0,1,0,1,1,0,1}, coord[2]={0.15,0.4},coord_m[2];
	int n_vertici=4;
	int puls_color[n_pulsanti]={rosso,giallo,blu,verde};
	char *menu[n_pulsanti]={"quit","cond.iniziale","Lungh.interaz.","probabilita"};
	colore=nero;
	for (int p=0;p<n_pulsanti;p++)
	{
		m_color(&puls_color[p]),m_polygon(vertici, &n_vertici);
		vertici[1]++,vertici[3]++,vertici[5]++,vertici[7]++;
		m_color(&colore),m_text(menu[p],coord),coord[1]++,m_flush();
	};
	//HO disegnato il menu
	int LI=5,A=0;
	string inizio="SSMMMRMMMMMRMRMMMMMMRMMMMMMRRRMMMMMRRMMMMRMRMMSSSSS";
	double seed=0;
	double PG=0.10,PC=0.25,PD=0.03,PN=0.15;
	automa_malattia_0(LI,inizio,PG,PC,PD,PN);
	while(A==0)
	{
		m_mouse(coord_m),m_which_window(&finestra);
		if (finestra==0)
		{
			if((coord_m[1]>3)&&(coord_m[1]<4))
			{
				int control=0;
				while (control==0)
				{
					cout<<"Immetti la probabilità di guarigione per un malato vicino a un medico: ",cin>>PG;
					if (PG<0||PG>1) {system("clear"),cerr<<"P deve stare tra 0 e 1!!"<<endl,control=0;}
					else control =1;
				};
				control=0;
				while (control==0)
				{
					cout<<"Immetti la probabilità di contagio per un sano vicino a un malato: ",cin>>PC;
					if (PC<0||PC>1) {system("clear"),cerr<<"P deve stare tra 0 e 1!!"<<endl,control=0;}
					else control =1;
				};
				control=0;
				while (control==0)
				{
					cout<<"Immetti la probabilità di morte per un malato (si sappia che malati vicino a morti muoiono piu facilmente in quanto si ritiene siano infetti di un ceppo letale!): ",cin>>PD;
					if (PD<0||PD>1) {system("clear"),cerr<<"P deve stare tra 0 e 1!!"<<endl,control=0;}
					else control =1;
				};
				control=0;
				while (control==0)
				{
					cout<<"Immetti la probabilità di natalità dove è morto qualcuno (rimpiazzamento popolazione): ",cin>>PN;
					if (PN<0||PN>1) {system("clear"),cerr<<"P deve stare tra 0 e 1!!"<<endl,control=0;}
					else control =1;
				};
				finestra=1,m_close(&finestra),automa_malattia_0(LI,inizio,PG,PC,PD,PN);
			};
			if((coord_m[1]>2)&&(coord_m[1]<3))
			{
				int control=0;
				while (control==0)
				{
					cout<<"Immetti la lunghezza d'interazione (intero positivo): ",cin>>LI;
					if (LI<1) {system("clear"),cerr<<"N deve valere almeno 1!!"<<endl,control=0;}
					else control =1;
				};
				finestra=1,m_close(&finestra),automa_malattia_0(LI,inizio,PG,PC,PD,PN);
			};
			if ((coord_m[1]>1)&&(coord_m[1]<2))
			{
				system("clear");
				cout<<"Legenda\nS-sano\nM-malato\nR-medico\nX-medico malato(non cura finchè lo è)\nD-deceduto\nI-immune\nImmetti la condizione iniziale: ", cin>>inizio,finestra=1,m_close(&finestra),automa_malattia_0(LI,inizio,PG,PC,PD,PN);
			}; //autogenerate
			if ((coord_m[1]<1)&&(coord_m[1]>0)) A=1;
		}
	};
	
	m_close(&(finestra=0)),m_close(&(finestra=1));
	return EXIT_SUCCESS;
}
int automa_malattia_0(int L,string a,float PG,float PC,float PD,float PN)
{
	ofstream fout ("automa.txt");
	int lunghezza=a.size(),guarigione,contagio,morire,k;//L lungh.d'interaz., varie probabilità che succeda un evento
	string b=a;
	double P,Pc,Pg,Pd;//PG probabilità di guarire con un medico vicino, PM probabilità di contagio con un malato vicino, PD probabilità di morte, di malattia, PN probabilità di nascita: occupazione del posto di un morto che è vacante.
	int  tempolimite=60;
	srand((time(0)));
	fout<<b<<endl;//stampa lo stato iniziale
	for (int istante=1; istante<=tempolimite;istante++)
	{
		for (int i=0; i<lunghezza;i++)
		{
			if (a[i]=='S')
			{
				guarigione=1,contagio=0,morire=0;
				for(int j=-L;j<=L;j++)
				{
					k=i+j;
					if(k<0)k=lunghezza+k;
					else if (k>lunghezza) k=k-lunghezza;
					if (a[k]=='M') contagio+=1;
					else if (a[k]=='R') {}
					else if (a[k]=='I') {}
					else if (a[k]=='D') {}
					else if (a[k]=='S') {}
					else if (a[k]=='X') contagio+=1;
				};
				Pc=(float)contagio*PC,Pg=(float)guarigione*PG,Pd=(float)morire*PD,P=Pc+Pg+Pd;
				float sorte=(float)rand()/RAND_MAX;
				if(P<1) P=1;
				if(0<=sorte&&sorte<=(Pc/P)) b[i]='M';
				else if((Pc/P)<=sorte&&sorte<=(Pc/P+Pd/P)) b[i]='D';
				else if((Pc/P+Pd/P)<=sorte&&sorte<=(Pc/P+Pd/P+Pg/P)) b[i]='S';
			}
			else if (a[i]=='M')
			{
				guarigione=1,contagio=0,morire=1;//di default un malato potrebbe morire!!! ma potrebbe anche curarsi
				for(int j=-L;j<=L;j++)
				{
					k=i+j;
					if(k<0)k=lunghezza+k;
					else if (k>lunghezza) k=k-lunghezza;
					if (a[i+j]=='M') {}
					else if (a[k]=='R') {guarigione+=1;}
					else if (a[k]=='I') {}
					else if (a[k]=='D') {morire+=1;}//ceppo letale,si trova in mezzo a dei morti!!!
					else if (a[k]=='S') {}
					else if (a[k]=='X') {}
				};
				Pc=(float)contagio*PC,Pg=(float)guarigione*PG,Pd=(float)morire*PD,P=Pc+Pg+Pd;
				float sorte=(float)rand()/RAND_MAX;
				if(P<1) P=1;
				if(0<=sorte&&sorte<=(Pc/P)) b[i]='M';
				else if((Pc/P)<=sorte&&sorte<=(Pc/P+Pd/P)) b[i]='D';
				else if((Pc/P+Pd/P)<=sorte&&sorte<=(Pc/P+Pd/P+Pg/P)) b[i]='S';
			}
			else if (a[i]=='R')
			{
				guarigione=0,contagio=0,morire=0;
				for(int j=-L;j<=L;j++)
				{
					k=i+j;
					if(k<0)k=lunghezza+k;
					else if (k>lunghezza) k=k-lunghezza;
					if (a[i+j]=='M') {contagio+=1;} 
					else if (a[k]=='R') {}
					else if (a[k]=='I') {}
					else if (a[k]=='D') {}
					else if (a[k]=='S') {}
					else if (a[k]=='X') {contagio+=1;}
				};
				Pc=(float)contagio*PC,Pg=(float)guarigione*PG,Pd=(float)morire*PD,P=Pc+Pg+Pd;
				float sorte=(float)rand()/RAND_MAX;
				if(P<1) P=1;
				if(0<=sorte&&sorte<=(Pc/P)) b[i]='X';
				else if((Pc/P)<=sorte&&sorte<=(Pc/P+Pd/P)) b[i]='D';
				else if((Pc/P+Pd/P)<=sorte&&sorte<=(Pc/P+Pd/P+Pg/P)) b[i]='R';
			}
			else if (a[i]=='X')
			{
				guarigione=1,contagio=0,morire=0;//di default un medico malato potrebbe morire!!!ma potrebbe anche curarsi
				for(int j=-L;j<=L;j++)
				{
					k=i+j;
					if(k<0)k=lunghezza+k;
					else if (k>lunghezza) k=k-lunghezza;
					if (a[i+j]=='M') {}
					else if (a[k]=='R') {guarigione+=1;}
					else if (a[k]=='I') {}
					else if (a[k]=='D') {morire+=1;}//ceppo letale,si trova in mezzo a dei morti!!!
					else if (a[k]=='S') {}
					else if (a[k]=='X') {guarigione+=1;}//un medico si puo curare da solo anche da malato
				};
				Pc=(float)contagio*PC,Pg=(float)guarigione*PG,Pd=(float)morire*PD,P=Pc+Pg+Pd;
				float sorte=(float)rand()/RAND_MAX;
				if(P<1) P=1;
				if(0<=sorte&&sorte<=(Pc/P)) b[i]='X';
				else if((Pc/P)<=sorte&&sorte<=(Pc/P+Pd/P)) b[i]='D';
				else if((Pc/P+Pd/P)<=sorte&&sorte<=(Pc/P+Pd/P+Pg/P)) b[i]='R';
			}
			else if (a[i]=='I'){}
			else if (a[i]=='D')
			{
				float sorte=(float)rand()/RAND_MAX;
				if(0<=sorte&&sorte<=PN) b[i]='S';
			};
		};
		a=b,fout<<b<<endl;
	};
	fout.close();
	plot_malattia(lunghezza,tempolimite,L,PG,PC,PD,PN);
	return 0;
}

int plot_malattia (int N,int iter,int lungint,float PG,float PC,float PD,float PN)
{
	system("clear");
	cout<<"probabilità di contagio per un sano vicino a un malato (per ogni malato) - "<<PC<<"\nprobabilità di guarigione per un malato vicino a un dottore (per ogni dottore) - "<<PG<<"\nprobabilità di morte per un malato (aumenta se ci sono morti in vicinanze) - "<<PD<<"\nprobabilità di natalità al posto di un morto - "<<PN<<"\nLunghezza automa - "<<N<<"\niterate - "<<iter<<"\nlunghezza d'interazione - "<<lungint<<"\n\nLegenda:\ngiallo-sano\nrosso-malato\nblu-medico\ncyan-medico malato(non cura finchè lo è)\nnero-deceduto\nmagenta-immune\n";
	char *t="Automa cellulare";
	ifstream file ("automa.txt");
	char automa;
	double xmin =-double(N),ymin=-double(iter+1),xmax=1.0,ymax=1.0, grandezza=1;
	int K=(100 << 24) + (100 << 16) + (100 << 8) + 100,s=3,L=(100<<24)+(N*14<<11)+14*(iter*+1);
	m_window (&K,&L);
	int finestra =2,colore=nero;
	m_use_as_pixmap (&finestra,&colore);
	m_frame(&xmin, &ymin, &xmax, &ymax);
	double coord[2];
	m_size_symbol(&grandezza);
	m_symbol(&s);
	for(int y=0;y>-(iter+1);y--)
	{for(int x=0;x>-N;x--)
	{
		coord[0]=x,coord[1]=y;
		file>>automa; 
		m_move(coord);
		if(automa=='R'){colore=blu,m_color(&colore),m_point(coord);}
		else if(automa=='D'){colore=nero,m_color(&colore),m_point(coord);}
		else if(automa=='I'){colore=magenta,m_color(&colore),m_point(coord);}
		else if(automa=='S'){colore=giallo,m_color(&colore),m_point(coord);}
		else if(automa=='M'){colore=rosso,m_color(&colore),m_point(coord);}
		else if(automa=='X'){colore=cyan,m_color(&colore),m_point(coord);}
		else;
	}
	}
	file.close();
	return EXIT_SUCCESS; 
}


