#include <FPT.h>
#include <D3d_matrix.h>
#include <stdlib.h>
#include <math.h>


typedef
struct {
  int objnum;
  int polynum;
  double dist;
}
  DISTANCES ;

DISTANCES painters[9000];
int numobjects;
int numpoints[10];
double x[10][9000], y[10][9000], z[10][9000];
int numpolys[10];
int psize[10][6000];
int con[10][6000][20];
double ambientConst, lightLoc[3], diffuseConst, intensity, hither=1;
double inherentRGB[3], actualRGB[3];

//==============================================================================

int  Clip_Polygon_Against_Plane(
				double a,double b,double c,double d, 
                  double *polyx, double *polyy, double *polyz, int size,
				double *resx, double *resy, double *resz)

// Clip polygon against the line ax + by + c = 0,
// where ax + by + c < 0 is considered IN.
// Incoming poly defined in arrays  polyx, polyy  with numverts = size.
// Clipped result values are stored in arrays  resx, resy,
// The numverts of the clipped result is returned as value of the function.

{
  int num,i,j ;
  double x1,y1,z1,x2,y2,z2,x21,y21,z21,den,t,xintsct,yintsct,zintsct ;
  double s1,s2;
  num = 0;

  for (i = 0 ; i < size ; i++) {
     j = (i + 1) % size ;

     // load up polygon to be clipped
     x1 = polyx[i] ; y1 = polyy[i] ; z1 = polyz[i] ;
     x2 = polyx[j] ; y2 = polyy[j] ; z2 = polyz[j] ;

     // clip line segment (x1,y1)-(x2,y2) against line
     s1 = (a*x1 + b*y1 + c*z1 + d) ;
     s2 = (a*x2 + b*y2 + c*z2 + d) ;

     if ((s1 >= 0) && (s2 >= 0)) {
        // out to out, do nothing
     } else if ((s1 < 0) && (s2 < 0)) {
        // in to in
       resx[num] = x2 ; resy[num] = y2 ; resz[num] = z2; num++; 
     } else {
        // one is in, the other out, so find the intersection

       x21 = x2 - x1 ; y21 = y2 - y1 ; z21 = z2 - z1;
       den = a*x21 + b*y21 + c*z21;
       if (den == 0) continue ; // do nothing-should never happen
       t = -(a*x1 + b*y1 + c*z1 + d)/den ;
       xintsct = x1 + t*x21 ;
       yintsct = y1 + t*y21 ;
       zintsct = z1 + t*z21;

       if (s1 < 0) { 
	 // in to out
	 resx[num] = xintsct ; resy[num] = yintsct ; resz[num] = zintsct; num++ ; 
       } else  {
	 // out to in
	 resx[num] = xintsct ; 
	 resy[num] = yintsct ; 
	 resz[num] = zintsct; 
	 num++ ;
	 resx[num] = x2 ; 
	 resy[num] = y2 ; 
	 resz[num] = z2; 
	 num++ ;
       }
       
     }
     
     
  } // end for i
  return num ;  // return size of the result poly
}

//==============================================================================

int  Clip_Polygon_Against_Pyramid (double *px,  double *py, double *pz, int psize){
  double nx[100],ny[100],nz[100], yonder=400 ;

   //top
   psize = Clip_Polygon_Against_Plane (0,1,-.57,0,
				      px,py,pz,psize,
				      nx,ny,nz) ;
   //bottom
   psize = Clip_Polygon_Against_Plane (0, -1,-.57,0,
				      nx,ny,nz,psize,
				      px,py,pz) ;
   //hither
   psize = Clip_Polygon_Against_Plane (0,0,-1,hither,
   				      px,py,pz,psize,
   				      nx,ny,nz) ;
   //yonder
   psize = Clip_Polygon_Against_Plane (0,0,1,-yonder,
   				      nx,ny,nz,psize,
  				      px,py,pz) ;
   //right
   psize = Clip_Polygon_Against_Plane (1,0,-.57,0,
				      px,py,pz,psize,
				      nx,ny,nz) ;
   //left
   psize = Clip_Polygon_Against_Plane (-1,0,-.57,0,
				      nx,ny,nz,psize,
				      px,py,pz) ;
   return psize ;
}
//=========================================================================
int readobject(FILE *f, int onum){
  int i, k, j;

  fscanf(f, "%d", &numpoints[onum]);
    
    for(k=0; k < numpoints[onum]; k++){
      fscanf(f, "%lf %lf %lf", &x[onum][k], &y[onum][k], &z[onum][k]);
    }
    fscanf(f, "%d", &numpolys[onum]);
    
    for(k=0; k < numpolys[onum]; k++){
      fscanf(f, "%d", &psize[onum][k]);
      
      for(j=0; j < psize[onum][k]; j++){
	fscanf(f, "%d", &con[onum][k][j]);
      }
    }   

}

//==============================================================================

int compare (const void *p, const void *q)
{
  DISTANCES *a, *b ;

  a = (DISTANCES*)p ;
  b = (DISTANCES*)q ;

  if  (((*a).dist) < ((*b).dist)) return -1 ;
  else if (((*a).dist) > ((*b).dist)) return 1 ;
  else return 0 ;
}

//==============================================================================

void lightModel(double polyX[], double polyY[], double polyZ[]){
  double unitLoc[3], overheadVector[3], vectorA[3], vectorB[3];
  double nlDotProduct, neDotproduct, erDotProduct, specular, diffuse, specPower=75;
  double lightDistance, originDistance, overHeadDistance;
  double LX, LY, LZ, NX, NY, NZ, RX, RY, RZ, EX, EY, EZ, checkValue, shading;
  int i ; 

  checkValue = ambientConst + diffuseConst ;
//calculating unit vector L
  LX = pow((lightLoc[0] - polyX[0]), 2) ;
  LY = pow((lightLoc[1] - polyY[0]), 2) ;
  LZ = pow((lightLoc[2] - polyZ[0]), 2) ;

  lightDistance = sqrt(LX + LY + LZ) ;

  LX = (lightLoc[0] - polyX[0])/lightDistance ;
  LY = (lightLoc[1] - polyY[0])/lightDistance ;
  LZ = (lightLoc[2] - polyZ[0])/lightDistance ;

//calculating vector of first line in poly
  vectorA[0] = polyX[1] - polyX[0] ;
  vectorA[1] = polyY[1] - polyY[0] ;
  vectorA[2] = polyZ[1] - polyZ[0] ;
//calculating vector of second line in poly
  vectorB[0] = polyX[2] - polyX[0] ;
  vectorB[1] = polyY[2] - polyY[0] ;
  vectorB[2] = polyZ[2] - polyZ[0] ;

//cross product of a and b to get overhead vector
  D3d_x_product(overheadVector, vectorA, vectorB) ;

//calculating unit vector N
  NX = pow((overheadVector[0] ), 2) ;
  NY = pow((overheadVector[1] ), 2) ;
  NZ = pow((overheadVector[2] ), 2) ;

  overHeadDistance = sqrt(NX + NY + NZ) ;

  NX = (overheadVector[0] )/overHeadDistance ;
  NY = (overheadVector[1] )/overHeadDistance ;
  NZ = (overheadVector[2] )/overHeadDistance ;

//finding dot product of N and L
  nlDotProduct = (NX*LX) + (NY*LY) + (NZ*LZ) ;


//Finding unit vector E
  EX = pow((polyX[0]), 2) ;
  EY = pow((polyY[0]), 2) ;
  EZ = pow((polyZ[0]), 2) ;

  originDistance = sqrt(EX + EY + EZ);

  EX = -polyX[0]/originDistance ;
  EY = -polyY[0]/originDistance ;
  EZ = -polyZ[0]/originDistance ;

//finding dot product of N and E

  neDotproduct = (NX*EX) + (NY*EY) + (NZ * EZ) ;


  if  (nlDotProduct * neDotproduct < 0){

    intensity = ambientConst ;
    //printf("ambient\n") ;

  } else if(nlDotProduct < 0 && neDotproduct < 0){
    NX = NX * -1 ;
    NY = NY * -1 ;
    NZ = NZ * -1 ;

    nlDotProduct *= -1 ;
    neDotproduct *= -1 ;
  }



//finding unit vector R
  RX = (2 * nlDotProduct * NX) - LX ; 
  RY = (2 * nlDotProduct * NY) - LY ;
  RZ = (2 * nlDotProduct * NZ) - LZ ;

//dot product of R and E gives angle beta

/*  double llen,nlen,elen,rlen ;

  llen = sqrt(LX*LX + LY*LY + LZ*LZ) ;
  nlen = sqrt(NX*NX + NY*NY + NZ*NZ) ;
  elen = sqrt(EX*EX + EY*EY + EZ*EZ) ;
  rlen = sqrt(RX*RX + RY*RY + RZ*RZ) ;

  printf("llen = %lf  nlen = %lf  elen = %lf   rlen = %lf\n",llen, nlen, elen,rlen) ;
  printf("nlDotProduct = %lf erDotProduct = %lf   ", nlDotProduct, erDotProduct) ;
  */

  erDotProduct = (RX*EX) + (RY*EY) + (RZ*EZ) ;

//diffuse
  diffuse = diffuseConst * nlDotProduct;

  specular = (1 - ambientConst - diffuseConst)* pow(erDotProduct, specPower) ;

  intensity = ambientConst + diffuse + specular ;

  if(intensity <= checkValue){
    shading = intensity/checkValue ;
    actualRGB[0] = shading * inherentRGB[0] ;
    actualRGB[1] = shading * inherentRGB[1] ;
    actualRGB[2] = shading * inherentRGB[2] ;
  }else{
    shading = (intensity - checkValue)/(1 - checkValue) ;
    actualRGB[0] = inherentRGB[0] + (shading * (1 - inherentRGB[0])) ;
    actualRGB[1] = inherentRGB[1] + (shading * (1 - inherentRGB[1])) ;
    actualRGB[2] = inherentRGB[2] + (shading * (1 - inherentRGB[2])) ;
  }
  return ;
}

//==============================================================================

void drawobject(){
  int i, k, m, f, j;
  double tmpX[9000], tmpY[9000], xx, yy, zz, halfAngle, lightX[100], lightY[100], lightZ[100];
  int onum,pnum,n,ntemp ;

  m = 0 ;
  for(i = 0; i < numobjects; i++){
    for(k = 0; k < numpolys[i]; k++){
      painters[m].objnum = i ;
      painters[m].polynum = k ;
      painters[m].dist = z[i][con[i][k][0]];
      m++ ;
    }
  }

  qsort (painters, m, sizeof(DISTANCES), compare) ;
  halfAngle = 30*M_PI/180 ;
  
  for(i = m-1; i >= 0; i--){

    onum = painters[i].objnum ;
    pnum = painters[i].polynum ;
    n =  psize[onum][pnum] ;

    for(k = 0; k < n ; k++){
      f = con[onum][pnum][k] ;
      lightX[k] = x[onum][f] ;
      lightY[k] = y[onum][f] ;
      lightZ[k] = z[onum][f] ;
    }


    n = Clip_Polygon_Against_Pyramid(lightX, lightY, lightZ, n);
    if (n == 0) continue ;


    for(k = 0; k < n ; k++){
      xx = lightX[k] ;
      yy = lightY[k] ;
      zz = lightZ[k] ;

      tmpX[k] = (300/(tan(halfAngle))) * (xx/zz) + 300 ;
      tmpY[k] = (300/(tan(halfAngle))) * (yy/zz) + 300 ;
    }

    lightModel(lightX, lightY, lightZ);
    G_rgb(actualRGB[0], actualRGB[1], actualRGB[2]) ;
    G_fill_polygon(tmpX, tmpY, n);
  }

}

//=============================================================================

int boundingbox(int counter){
  int i=0;
  double lowx , highx, lowy, highy, lowz, highz, sfx, sfy, sfz, sf ;
  double  cx, cy, cz, m[4][4], mi[4][4];
  lowx = x[counter][i];
  lowy = y[counter][i];
  lowz = z[counter][i];
  highx = x[counter][i];
  highy = y[counter][i];
  highz = z[counter][i];

    for(i=0; i<numpoints[counter]; i++){

      if(x[counter][i]<lowx){
	lowx = x[counter][i];
      }
      if(x[counter][i]>highx){
	highx = x[counter][i];
      }
      if(y[counter][i]<lowy){
	lowy = y[counter][i];
      }
      if(y[counter][i]>highy){
	highy = y[counter][i];
      }
      if(z[counter][i]<lowz){
	lowz = z[counter][i];
      }
      if(z[counter][i]>highz){
	highz = z[counter][i];
      }
      
    }

  cx = (lowx + highx)/2 ;
  cy = (lowy + highy)/2 ;
  cz = (lowz + highz)/2 ;

  D3d_make_identity(m) ;
  D3d_make_identity(mi) ;
  D3d_translate(m, mi, -cx, -cy, -cz) ;

  sfx = 1/(highx - lowx) ;
  sfy = 1/(highy - lowy) ;
  sfz = 1/(highz - lowz) ;
  if (sfy > sfx && sfz > sfx) {
    sf = sfx ;
  } else if(sfx > sfy && sfz > sfy){
    sf = sfy;
  }else{
    sf = sfz;
  }

  D3d_scale(m, mi, sf, sf, sf) ;
  // printf("%d\n", counter);
  D3d_mat_mult_points (x[counter],y[counter], z[counter],  
		       m, x[counter], y[counter], z[counter], numpoints[counter]) ;
  // printf("through bounding box\n");
  // adding the centerpoints
  x[counter][numpoints[counter]] = 0;
  y[counter][numpoints[counter]] = 0;
  z[counter][numpoints[counter]] = 0;
  return 1;

}

//==============================================================================

int main(int argc, char **argv){

  FILE *f;
  int i, counter, objectSelected;
  double  m[4][4], mi[4][4], tmpx[100], tmpy[100];

  if (argc < 2){
      printf("usage: program file name");
  }

  numobjects = argc - 1 ;

  for(i=0; i < numobjects; i++){
    
    f = fopen(argv[i+1], "r");
    
    if(f == NULL){
      printf("Can't open %s\n", argv[i]);
      exit(0);
    }
    readobject(f,i);
    boundingbox(i);
  }
 
  //takes in all light info
  printf("Please input x, y and z coordinates for the light\n");
  scanf("%lf %lf %lf", &lightLoc[0], &lightLoc[1], &lightLoc[2]);
  printf("Please input the ambient light constant\n");
  scanf("%lf", &ambientConst);
  printf("Please input the diffuse constant\n");
  scanf("%lf", &diffuseConst);
  printf("please input the inherent color of each object\n");
  scanf("%lf %lf %lf", &inherentRGB[0], &inherentRGB[1], &inherentRGB[2]);
  
  G_init_graphics(600, 600);
  G_rgb(0,0,0);
  G_clear();
  drawobject();
  //  printf("%d\n", numpoints[counter]);

  counter = 0 ;
  while(counter != 'q'){
    counter = G_wait_key();
    G_rgb(0,0,0);
    switch(counter){

      //move away in z axis
    case 'a':
      D3d_make_identity(m);
      D3d_make_identity(mi);
      D3d_translate(m, mi, 0, 0, .05);
      D3d_mat_mult_points(x[objectSelected], y[objectSelected], 
			  z[objectSelected], m, x[objectSelected], 
			  y[objectSelected], z[objectSelected], 
			  numpoints[objectSelected]+1);

      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //move closer in z axis
    case 'z':
      D3d_make_identity(m);
      D3d_make_identity(mi);
      D3d_translate(m, mi, 0, 0, -.05);
      D3d_mat_mult_points(x[objectSelected], y[objectSelected], 
			  z[objectSelected], m, x[objectSelected], 
			  y[objectSelected], z[objectSelected], 
			  numpoints[objectSelected]+1);
      
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //move up in y
    case 's':
      D3d_make_identity(m);
      D3d_make_identity(mi);
      D3d_translate(m, mi, 0, .05, 0);
      D3d_mat_mult_points(x[objectSelected], y[objectSelected], 
			  z[objectSelected], m, x[objectSelected], 
			  y[objectSelected], z[objectSelected], 
			  numpoints[objectSelected]+1);
      
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //move down in y
    case 'x':
      D3d_make_identity(m);
      D3d_make_identity(mi);
      D3d_translate(m, mi, 0, -.05, 0);
      D3d_mat_mult_points(x[objectSelected], y[objectSelected], 
			  z[objectSelected], m, x[objectSelected], 
			  y[objectSelected], z[objectSelected], 
			  numpoints[objectSelected]+1);
      
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //move to the right
    case 'd':
      D3d_make_identity(m);
      D3d_make_identity(mi);
      D3d_translate(m, mi, .05, 0, 0);
      D3d_mat_mult_points(x[objectSelected], y[objectSelected], 
			  z[objectSelected], m, x[objectSelected], 
			  y[objectSelected], z[objectSelected], 
			  numpoints[objectSelected]+1);
      
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //move to the left
    case 'c':
      D3d_make_identity(m);
      D3d_make_identity(mi);
      D3d_translate(m, mi, -.05, 0, 0);
      D3d_mat_mult_points(x[objectSelected], y[objectSelected], 
			  z[objectSelected], m, x[objectSelected], 
			  y[objectSelected], z[objectSelected], 
			  numpoints[objectSelected]+1);
      
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //right rotation along z axis
    case 'f':
      D3d_make_identity(m);
      D3d_make_identity(mi);
      D3d_translate(m, mi, -x[objectSelected][numpoints[objectSelected]], 
		    -y[objectSelected][numpoints[objectSelected]], 
		    -z[objectSelected][numpoints[objectSelected]]);
      D3d_rotate_z(m, mi, .05);
      D3d_translate(m, mi, x[objectSelected][numpoints[objectSelected]], 
		    y[objectSelected][numpoints[objectSelected]],
		    z[objectSelected][numpoints[objectSelected]]);
      D3d_mat_mult_points(x[objectSelected], y[objectSelected], 
			  z[objectSelected], m, x[objectSelected], 
			  y[objectSelected], z[objectSelected], 
			  numpoints[objectSelected]+1);
      
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;
      
      //left rotation along z axis
    case 'v':
      D3d_make_identity(m);
      D3d_make_identity(mi);
      D3d_translate(m, mi, -x[objectSelected][numpoints[objectSelected]], 
		    -y[objectSelected][numpoints[objectSelected]], 
		    -z[objectSelected][numpoints[objectSelected]]);
      D3d_rotate_z(m, mi, -.05);
      D3d_translate(m, mi, x[objectSelected][numpoints[objectSelected]], 
	      y[objectSelected][numpoints[objectSelected]],
	      z[objectSelected][numpoints[objectSelected]]);
      D3d_mat_mult_points(x[objectSelected], y[objectSelected], 
			  z[objectSelected], m, x[objectSelected], 
			  y[objectSelected], z[objectSelected], 
			  numpoints[objectSelected]+1);
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //rotation along y axis to the left
    case 'g':
      D3d_make_identity(m);
      D3d_make_identity(mi);
      D3d_translate(m, mi, -x[objectSelected][numpoints[objectSelected]], 
		    -y[objectSelected][numpoints[objectSelected]], 
		    -z[objectSelected][numpoints[objectSelected]]);
      D3d_rotate_y(m, mi, .05);
      D3d_translate(m, mi, x[objectSelected][numpoints[objectSelected]], 
		    y[objectSelected][numpoints[objectSelected]],
		    z[objectSelected][numpoints[objectSelected]]);
      D3d_mat_mult_points(x[objectSelected], y[objectSelected], 
			  z[objectSelected], m, x[objectSelected], 
			  y[objectSelected], z[objectSelected], 
			  numpoints[objectSelected]+1);
      
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //rotation along y axis to the right
    case 'b':
      D3d_make_identity(m);
      D3d_make_identity(mi);
      D3d_translate(m, mi, -x[objectSelected][numpoints[objectSelected]], 
		    -y[objectSelected][numpoints[objectSelected]], 
		    -z[objectSelected][numpoints[objectSelected]]);
      D3d_rotate_y(m, mi, -.05);
      D3d_translate(m, mi, x[objectSelected][numpoints[objectSelected]], 
	      y[objectSelected][numpoints[objectSelected]],
	      z[objectSelected][numpoints[objectSelected]]);
      D3d_mat_mult_points(x[objectSelected], y[objectSelected], 
			  z[objectSelected], m, x[objectSelected], 
			  y[objectSelected], z[objectSelected], 
			  numpoints[objectSelected]+1);
      
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //rotation along x axis upwards
    case 'h':
      D3d_make_identity(m);
      D3d_make_identity(mi);
      D3d_translate(m, mi, -x[objectSelected][numpoints[objectSelected]], 
		    -y[objectSelected][numpoints[objectSelected]], 
		    -z[objectSelected][numpoints[objectSelected]]);
      D3d_rotate_x(m, mi, .05);
      D3d_translate(m, mi, x[objectSelected][numpoints[objectSelected]], 
	      y[objectSelected][numpoints[objectSelected]],
	      z[objectSelected][numpoints[objectSelected]]);
      D3d_mat_mult_points(x[objectSelected], y[objectSelected], 
			  z[objectSelected], m, x[objectSelected], 
			  y[objectSelected], z[objectSelected], 
			  numpoints[objectSelected]+1);
      
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //rotation along x axis downwards
    case 'n':
      D3d_make_identity(m);
      D3d_make_identity(mi);
      D3d_translate(m, mi, -x[objectSelected][numpoints[objectSelected]], 
		    -y[objectSelected][numpoints[objectSelected]], 
		    -z[objectSelected][numpoints[objectSelected]]);
      D3d_rotate_x(m, mi, -.05);
      D3d_translate(m, mi, x[objectSelected][numpoints[objectSelected]], 
		    y[objectSelected][numpoints[objectSelected]],
		    z[objectSelected][numpoints[objectSelected]]);
      D3d_mat_mult_points(x[objectSelected], y[objectSelected], 
			  z[objectSelected], m, x[objectSelected], 
			  y[objectSelected], z[objectSelected], 
			  numpoints[objectSelected]+1);
      
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //these case statements move the light
      //moves the light right
    case 'u':
      lightLoc[0] = lightLoc[0] + 10;
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      // moves the light up
    case 'i':
      lightLoc[1] = lightLoc[1] + 10;
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //moves the light forward
    case 'o':
      lightLoc[2] = lightLoc[2] + 10;
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //moves the light left
    case 'j':
      lightLoc[0] = lightLoc[0] - 10;
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //moves the light down
    case 'k':
      lightLoc[1] = lightLoc[1] - 10;
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //moves the light back
    case 'l':
      lightLoc[2] = lightLoc[2] - 10;
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;
     
      //moves hither forward
    case 't':
      hither = hither+.1;
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //moves hither back
    case'y':
      hither = hither-.1;
      G_rgb(0,0,0);
      G_clear();
      drawobject();
      break;

      //switches object selected to manipulate
    case '0':
      objectSelected = 0;
      break;
    case '1':
      objectSelected = 1;
      break;
    case '2':
      objectSelected = 2;
      break;
    case '3':
      objectSelected = 3;
      break;
    case '4':
      objectSelected = 4;
      break;
    case '5':
      objectSelected = 5;
      break;
    case '6':
      objectSelected = 6;
      break;
    case '7':
      objectSelected = 7;
      break;
    case '8':
      objectSelected = 8;
      break;
    case '9':
      objectSelected = 9;
      break;

      //checks to make sure key pressed is defined
    default:
      if(counter != 'q'){
      printf("you pressed a key that isn't defined, try again\n");
      }
      break;
    }
  }
}
