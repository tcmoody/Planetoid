#include <FPT.h>
#include <math.h>
#include <D3d_matrix.h>

double manMatrix[20][4][4], invManMatrix[20][4][4];
double reflectivity[20];
double inherentColor[20][3];
double light_in_eye_space[3] ;
double refractivity[20];
double refractiveIndex[20];
double backgroundRefractiveIndex=1;
int numObjects;

double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 50 ;

int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3])
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// s = location of start of ray (probably the eye)
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]

// return 1 if successful, 0 if error
{

  double len ;
  double N[3] ; 
  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
  if (len == 0) return 0 ;
  N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;

  double E[3] ;
  E[0] = s[0] - p[0] ; 
  E[1] = s[1] - p[1] ; 
  E[2] = s[2] - p[2] ; 
  len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
  if (len == 0) return 0 ;
  E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
  double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;

  double L[3] ;
  L[0] = light_in_eye_space[0] - p[0] ; 
  L[1] = light_in_eye_space[1] - p[1] ; 
  L[2] = light_in_eye_space[2] - p[2] ; 
  len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
  if (len == 0) return 0 ;
  L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;





  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
     // this needs to occur BEFORE you possibly jump to LLL below




  double intensity ;
  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = AMBIENT ; 
    goto LLL ;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // eye and light on same side but normal pointing "wrong" way
    N[0] *= (-1.0) ;    N[1] *= (-1.0) ;    N[2] *= (-1.0) ; 
    NdotL *= (-1.0) ;
    NdotE *= (-1.0) ;   // don't use NdotE below, probably should eliminate this
  }


  // ignore Blinn's variant
  double R[3] ; // Reflection vector of incoming light
  R[0] = 2*NdotL*N[0] - L[0] ;
  R[1] = 2*NdotL*N[1] - L[1] ;
  R[2] = 2*NdotL*N[2] - L[2] ;

  double EdotR = E[0]*R[0] + E[1]*R[1] + E[2]*R[2] ;

  double diffuse ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = MAX_DIFFUSE*NdotL ; }

  double specular ;
  if (EdotR <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPECPOW) ;}

  // printf("%lf %lf\n",diffuse,specular) ;
  intensity = AMBIENT + diffuse + specular ;



 LLL : ;

  double f,g ;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse ;
    argb[0] = f * irgb[0] ;
    argb[1] = f * irgb[1] ;
    argb[2] = f * irgb[2] ;
  } else {
    f = (intensity - max_ambient_and_diffuse) / 
                           (1.0 - max_ambient_and_diffuse) ;
    g = 1.0 - f ;
    argb[0] = g * irgb[0] + f ;
    argb[1] = g * irgb[1] + f ;
    argb[2] = g * irgb[2] + f ;
  }

  return 1 ;
}

//==========
//need to apply the transformations to the inverse slopes
void makeManMatrix(int count, double sx, double sy, double sz, double rx, double ry, double rz, double tx, double ty, double tz){
  int num=0;
  int tlist[9];
  double plist[9];

  tlist[num] = SX ; plist[num] =  sx ; num++ ;
  tlist[num] = SY ; plist[num] =  sy ; num++ ;
  tlist[num] = SZ ; plist[num] =  sz ; num++ ;
  tlist[num] = RX ; plist[num] =  rx ; num++ ;
  tlist[num] = RY ; plist[num] =  ry ; num++ ;
  tlist[num] = RZ ; plist[num] =  rz ; num++ ;
  tlist[num] = TX ; plist[num] =  tx ; num++ ;
  tlist[num] = TY ; plist[num] =  ty ; num++ ;
  tlist[num] = TZ ; plist[num] =  tz ; num++ ;


  D3d_make_movement_sequence_matrix(manMatrix[count], invManMatrix[count], num, tlist, plist);
}

int quadratic(double a, double b, double c, double t[2]){
	double under=(b*b)-(4*a*c);
        if (under < 0) return -1 ;
        under = sqrt(under) ;
	t[0] = (-b+under)/(2*a);
	t[1] = (-b-under)/(2*a);
	return 2 ;
}

void findReflectsAndRefracts(int numReflections, double s[3], double e[3], double color[3], int prevObj, int refractCounter){

  if(numReflections==0){ return;	}
  
  double pointOne[3], pointTwo[3] ;
  double tempPointOne[3], tempPointTwo[3] ;
  pointOne[0]=s[0]; pointOne[1]=s[1]; pointOne[2]=s[2] ;
  pointTwo[0]=e[0]; pointTwo[1]=e[1]; pointTwo[2]=e[2] ;
  
  double t[2] ;
  int i,sig ;
  double tsave[100] ;
  int osave[100] ;
  int numsaved = 0 ;
  for(i=0; i<numObjects; i++){ // make a var for num objects
    
    D3d_mat_mult_pt(tempPointOne, invManMatrix[i], pointOne);
    D3d_mat_mult_pt(tempPointTwo, invManMatrix[i], pointTwo);
    
    double a, b, c, p, q, r;
    p=tempPointTwo[0]-tempPointOne[0];
    q=tempPointTwo[1]-tempPointOne[1];
    r=tempPointTwo[2]-tempPointOne[2];
    
    a=(p*p)+(q*q)+(r*r);
    b=(2*tempPointOne[0]*p)+(2*tempPointOne[1]*q)+(2*tempPointOne[2]*r);
    c=(tempPointOne[0]*tempPointOne[0])+(tempPointOne[1]*tempPointOne[1])+(tempPointOne[2]*tempPointOne[2])-1;
    
    sig = quadratic(a, b, c, t);
    if (sig == -1) { continue ; }
    if (t[0] > 0) { tsave[numsaved] = t[0] ; osave[numsaved] = i ; numsaved++ ; }
    if (t[1] > 0) { tsave[numsaved] = t[1] ; osave[numsaved] = i ; numsaved++ ; }
    
  }
  
  if (numsaved == 0) { 
    color[0]=0; 
    color[1]=0; 
    color[2]=0; 
    return;
  }
  
  double minposT = 1e50 ;
  int minobj = -1 ;
  for (i = 0 ; i < numsaved ; i++) {
      if (tsave[i] < minposT) { minposT = tsave[i] ; minobj = osave[i] ; }
  }
  
  double intersectPointEyeSpace[3] ;
  double intersectPointObjSpace[3] ;
  
  intersectPointEyeSpace[0]=pointOne[0]+minposT*(pointTwo[0] - pointOne[0]) ;
  intersectPointEyeSpace[1]=pointOne[1]+minposT*(pointTwo[1] - pointOne[1]) ;
  intersectPointEyeSpace[2]=pointOne[2]+minposT*(pointTwo[2] - pointOne[2]) ;
  
  D3d_mat_mult_pt(intersectPointObjSpace, invManMatrix[minobj], intersectPointEyeSpace);
  
  double normalObjSpace[3] ;
  normalObjSpace[0] = 2*intersectPointObjSpace[0] ;
  normalObjSpace[1] = 2*intersectPointObjSpace[1] ;
  normalObjSpace[2] = 2*intersectPointObjSpace[2] ;
  
  //FIX THIS
  double normalEyeSpace[3] ;
  normalEyeSpace[0] = invManMatrix[minobj][0][0]*normalObjSpace[0] + invManMatrix[minobj][1][0]*normalObjSpace[1] + invManMatrix[minobj][2][0]*normalObjSpace[2] ;
  normalEyeSpace[1] = invManMatrix[minobj][0][1]*normalObjSpace[0] + invManMatrix[minobj][1][1]*normalObjSpace[1] + invManMatrix[minobj][2][1]*normalObjSpace[2] ;
  normalEyeSpace[2] = invManMatrix[minobj][0][2]*normalObjSpace[0] + invManMatrix[minobj][1][2]*normalObjSpace[1] + invManMatrix[minobj][2][2]*normalObjSpace[2] ;  

  double len ;
  double normal[3] ;
  normal[0] = normalEyeSpace[0] ;
  normal[1] = normalEyeSpace[1] ;
  normal[2] = normalEyeSpace[2] ;
  len = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]) ;
  normal[0] /= len ;
  normal[1] /= len ;
  normal[2] /= len ;

   double nustart[3], nuend[3];
if(refractivity[minobj]!=0){
  double n1, n2;
  if(refractCounter==0 || (fmod(refractCounter, 2)==0)){
    n1=backgroundRefractiveIndex;
    n2=refractiveIndex[minobj];
  }else{
    n1=refractiveIndex[minobj];
    n2=backgroundRefractiveIndex;
  }
        
  double refractUnitVector[3], intersectUnitVector[3], nDotI=0;
  
  //getting intersect unit vector
  len=sqrt(pow(s[0]-e[0],2) + pow(s[1]-e[1],2) + pow(s[2]-e[2],2));
  intersectUnitVector[0]=(e[0]-s[0])/len;
  intersectUnitVector[1]=(e[1]-s[1])/len;
  intersectUnitVector[2]=(e[2]-s[2])/len;
  

  //dot product of intersect unit vector and normal unit vector
  double rnormal[3] ;
  rnormal[0] = normal[0] ;
  rnormal[1] = normal[1] ;
  rnormal[2] = normal[2] ;
  nDotI=(intersectUnitVector[0]*rnormal[0])+
        (intersectUnitVector[1]*rnormal[1])+
    (intersectUnitVector[2]*rnormal[2]);
  
  if(nDotI<0){
    nDotI = - nDotI ;
    rnormal[0] = -rnormal[0] ;
    rnormal[1] = -rnormal[1] ;
    rnormal[2] = -rnormal[2] ;
  }
  
  //calculating refraction unit vector

  double under=(1- pow(n1/n2,2)*(1-pow(nDotI,2)));
  if (under > 0) {
    under = sqrt(under) ;
    refractUnitVector[0]=((n1/n2)*intersectUnitVector[0])+((under-(n1/n2)*(nDotI))*rnormal[0]);
    refractUnitVector[1]=((n1/n2)*intersectUnitVector[1])+((under-(n1/n2)*(nDotI))*rnormal[1]);
    refractUnitVector[2]=((n1/n2)*intersectUnitVector[2])+((under-(n1/n2)*(nDotI))*rnormal[2]);

    nustart[0] = intersectPointEyeSpace[0] + 0.001*refractUnitVector[0] ;
    nustart[1] = intersectPointEyeSpace[1] + 0.001*refractUnitVector[1] ;
    nustart[2] = intersectPointEyeSpace[2] + 0.001*refractUnitVector[2] ;

    nuend[0] = intersectPointEyeSpace[0] + 0.002*refractUnitVector[0] ;
    nuend[1] = intersectPointEyeSpace[1] + 0.002*refractUnitVector[1] ;
    nuend[2] = intersectPointEyeSpace[2] + 0.002*refractUnitVector[2] ;
    findReflectsAndRefracts(numReflections-1, nustart, nuend, color, minobj, refractCounter+1);

    double tempColor[3];
    tempColor[0]=color[0]; tempColor[1]=color[1]; tempColor[2]=color[2];
  // irgb == inherent color of object (input to this function)
  // s = location of start of ray (probably the eye)
  // p = point on object (input to this function)
  // n = normal to the object at p (input to this function)
  // argb == actual color of object (output of this function)
  //call to lightmodel  
    Light_Model(inherentColor[minobj], s, intersectPointEyeSpace, normalEyeSpace, color);
    if(prevObj!=100){
      color[0]=tempColor[0]*refractivity[prevObj]+color[0]*refractivity[minobj];
      color[1]=tempColor[1]*refractivity[prevObj]+color[1]*refractivity[minobj];
      color[2]=tempColor[2]*refractivity[prevObj]+color[2]*refractivity[minobj];
    }else{
      color[0]=tempColor[0]*refractivity[minobj];
      color[1]=tempColor[1]*refractivity[minobj];
      color[2]=tempColor[2]*refractivity[minobj];
    }
  }
}else{

  double lvector[3] ;
  lvector[0] = pointOne[0] - pointTwo[0] ;
  lvector[1] = pointOne[1] - pointTwo[1] ;
  lvector[2] = pointOne[2] - pointTwo[2] ;
  len = sqrt(lvector[0]*lvector[0] + lvector[1]*lvector[1] + lvector[2]*lvector[2]) ;
  lvector[0] /= len ;
  lvector[1] /= len ;
  lvector[2] /= len ;
  
  double ndotl = normal[0]*lvector[0] + normal[1]*lvector[1] + normal[2]*lvector[2]  ;
  double reflection[3] ;
  reflection[0] = 2*ndotl*normal[0] - lvector[0] ;
  reflection[1] = 2*ndotl*normal[1] - lvector[1] ;
  reflection[2] = 2*ndotl*normal[2] - lvector[2] ;
  
  nustart[0] = intersectPointEyeSpace[0] + 0.01*reflection[0] ;
  nustart[1] = intersectPointEyeSpace[1] + 0.01*reflection[1] ;
  nustart[2] = intersectPointEyeSpace[2] + 0.01*reflection[2] ;
  
  nuend[0] = intersectPointEyeSpace[0] + 1.00*reflection[0] ;
  nuend[1] = intersectPointEyeSpace[1] + 1.00*reflection[1] ;
  nuend[2] = intersectPointEyeSpace[2] + 1.00*reflection[2] ;

  findReflectsAndRefracts(numReflections-1, nustart, nuend, color, minobj, refractCounter) ;

  double tempColor[3];
  tempColor[0]=color[0]; tempColor[1]=color[1]; tempColor[2]=color[2];
  // irgb == inherent color of object (input to this function)
  // s = location of start of ray (probably the eye)
  // p = point on object (input to this function)
  // n = normal to the object at p (input to this function)
  // argb == actual color of object (output of this function)
  //call to lightmodel  
  Light_Model(inherentColor[minobj], s, intersectPointEyeSpace, normalEyeSpace, color);

  if(prevObj!=100){
    color[0]+=tempColor[0]*reflectivity[prevObj] + color[0]*reflectivity[minobj];
    color[1]+=tempColor[1]*reflectivity[prevObj] + color[1]*reflectivity[minobj];
    color[2]+=tempColor[2]*reflectivity[prevObj] + color[2]*reflectivity[minobj];
  }else{
    color[0]+=tempColor[0]*reflectivity[minobj];
    color[1]+=tempColor[1]*reflectivity[minobj];
    color[2]+=tempColor[2]*reflectivity[minobj];
  }
}
return;
}

void drawAllObjects(){
  //ray starts from origin to pixel
  int i, j;
  double tan_degrees_of_half_angle = tan((30*M_PI)/180);
  double s[3], e[3], color[3], result=0;

  s[0]=0; s[1]=0; s[2]=0;
  e[2]= 300/tan_degrees_of_half_angle;

  for(i=-300; i< 300; i++){
    for (j=-300; j < 300; j++){
      e[0]=i ;
      e[1]=j ;
      findReflectsAndRefracts(3, s, e,  color, 100, 0);
      G_rgb(color[0], color[1], color[2]);
      G_point(i+300,j+300) ;
    }
  }
}

int main(){
  int key;
  G_init_graphics(600,600);
  double i=-200;
  int j=0; 
  int k=0;
  light_in_eye_space[0]=100; light_in_eye_space[1]=200; light_in_eye_space[2]= -100;
  numObjects=16;
  char prefix[100], filename[100];
  strncpy(prefix, "planetoid", 100);

  while(j<350){
    printf("j: %d\n", j);
    //lens
    makeManMatrix(0, 100, 100, 5, 0, 0, 0, i, 0, 200);
    //planet
    makeManMatrix(1, 400, 400, 400, 0, 0, 0, 0, 0, 15000);
    //moon
    makeManMatrix(2, 150, 150, 150, 0, 0, 0, 500, 500, 15000);
    //asteroids
    makeManMatrix(3, 45, 45, 45, 0, 0, 0, 290, -340, 14900);
    makeManMatrix(4, 50, 50, 50, 0, 0, 0, 450, -450, 15000);


    reflectivity[0]=0; 
    reflectivity[1]=.5; 
    reflectivity[2]=.5; 
    reflectivity[3]=.5; 
    reflectivity[4]=.5;
    reflectivity[5]=.5;


    refractivity[0]=.8; 
    refractivity[1]=0; 
    refractivity[2]=0; 
    refractivity[3]=0; 
    refractivity[4]=0;
    refractivity[5]=0;


    refractiveIndex[0]=4; 
    refractiveIndex[1]=1.4; 
    refractiveIndex[2]=1.4; 
    refractiveIndex[3]=1.4; 
    refractiveIndex[4]=1.4;
    refractiveIndex[5]=1.4;

        //lens
    inherentColor[0][0]=.7; inherentColor[0][1]=.7; inherentColor[0][2]=.7;
    //planet
    inherentColor[1][0]=0; inherentColor[1][1]=.3; inherentColor[2][2]=.4;
    //moon
    inherentColor[2][0]=.2; inherentColor[2][1]=.2; inherentColor[2][2]=.2;
    //asteroids
    inherentColor[3][0]=.2; inherentColor[3][1]=0; inherentColor[3][2]=0;
    inherentColor[4][0]=0; inherentColor[4][1]=0; inherentColor[4][2]=.2;


    //stars
    makeManMatrix(5, 800, 800, 800, 0, 0, 0, 30000, -30000, 400000);
    makeManMatrix(6, 800, 800, 800, 0, 0, 0, 30000, 30000, 400000);
    makeManMatrix(7, 800, 800, 800, 0, 0, 0, -30000, 30000, 400000);
    makeManMatrix(8, 800, 800, 800, 0, 0, 0, -30000, -30000, 400000);
    makeManMatrix(9, 800, 800, 800, 0, 0, 0, -50000, -10000, 400000);
    makeManMatrix(10, 800, 800, 800, 0, 0, 0, 65000, -30000, 400000);
    makeManMatrix(11, 800, 800, 800, 0, 0, 0, 35000, 65000, 400000);
    makeManMatrix(12, 800, 800, 800, 0, 0, 0, 15000, -30000, 400000);
    makeManMatrix(13, 800, 800, 800, 0, 0, 0, -5000, 61000, 400000);
    makeManMatrix(14, 800, 800, 800, 0, 0, 0, 40000, 0, 400000);
    makeManMatrix(15, 800, 800, 800, 0, 0, 0, 0, -45000, 400000);

    for(k=5; k<16; k++){
      inherentColor[k][0]=1; inherentColor[k][1]=1; inherentColor[k][2]=1;
      reflectivity[k]=.5;
      refractivity[k]=0;
      refractiveIndex[k]=1.4;
    }

    if(i>=-26){
      i+=.5;
    }else{
      i++;
    }
    drawAllObjects();
    sprintf(filename, "%s%04d", prefix, j);
    filename[13] = '.';
    filename[14] = 'x';
    filename[15] = 'w';
    filename[16] = 'd';
    G_save_image_to_file(filename);
    j++;
    G_rgb(1,1,1);
    G_clear();  
  }
}
