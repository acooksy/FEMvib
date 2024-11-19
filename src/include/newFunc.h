//
//  newFunc.h
//  GMAT
//
//  Created by PETER ZAJAC on 10/5/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

struct functionRot
{
    
    MatDoub xref;
    MatDoub xgeom;
    Int natom; 
    VecDoub mass;

    functionRot(MatDoub xre, MatDoub xg, Int na, VecDoub ma):xref(xre), xgeom(xg), natom(na), mass(ma) { }
    // note ee[0] ee[1] ee[2] hold the vector of rotation
    // and ee[3] is the angle 0 to 2*pi
    
    Doub operator()(VecDoub pokus)
    {
        
        
    int i,j,k;
    double sizeE=0.0;
    double pi = 3.1415926535897932384626433;
    double cs,sn,cc,sp,result=0.0;
    MatDoub rotstruct(natom,3);
    MatDoub angmom(natom,3);
    VecDoub vp(3);
    VecDoub tempdiff(3);
    VecDoub vecprod(3);
    VecDoub sum(3);
    
    
    // normalization
    for (i=0; i<3; i++) 
    {
        sizeE = sizeE  + pokus[i]*pokus[i];
    }
    sizeE = sqrt(sizeE);
    
    for (i=0; i<3; i++) 
    {
        pokus[i] = pokus[i]/sizeE;
    }    
    
    // setting ee[3] to be within 0 and pi
  //  while (abs(pokus[3]) > pi) 
    while (abs(pokus[3]) > pi)
    {
        if (pokus[3] < 0.0) 
        {
            pokus[3] = pokus[3] + pi;
        }
        if (pokus[3] > pi)
        {
            pokus[3] = pokus[3] - pi;
        }
    }
        
        if (pokus[3] < 0.0) 
        {
            pokus[3] = pokus[3] + pi;
        }
    
    // rotating the structure
    cs = cos(pokus[3]);
    sn = sin(pokus[3]);

    for (i=0; i<natom; i++ ) 
    {
        sp = 0.0;
        for (j=0; j<3; j++) 
        {
            sp = sp + xgeom[i][j]*pokus[j];
        }
        cc = (1.0-cs)*sp;
        vp[0] = pokus[1]*xgeom[i][2] - pokus[2]*xgeom[i][1];
        vp[1] = pokus[2]*xgeom[i][0] - pokus[0]*xgeom[i][2];
        vp[2] = pokus[0]*xgeom[i][1] - pokus[1]*xgeom[i][0];
        rotstruct[i][0] = cs*xgeom[i][0] + cc*pokus[0] + sn*vp[0];
        rotstruct[i][1] = cs*xgeom[i][1] + cc*pokus[1] + sn*vp[1];
        rotstruct[i][2] = cs*xgeom[i][2] + cc*pokus[2] + sn*vp[2];
    }
    // rotation is finished 
    
    // calculate angular momentum
/*    for (i=0; i<natom; i++) 
    {
    //   tempdiff[0] = rotstruct[i][0] - xref[i][0];
    //    tempdiff[1] = rotstruct[i][1] - xref[i][1];
    //    tempdiff[2] = rotstruct[i][2] - xref[i][2];
    
        
     // vector product between xrot and tempdiff
      //  vecprod[0] = rotstruct[i][1]*tempdiff[2] - rotstruct[i][2]*tempdiff[1];
      //  vecprod[1] = rotstruct[i][2]*tempdiff[0] - rotstruct[i][0]*tempdiff[2];
      //  vecprod[2] = rotstruct[i][0]*tempdiff[1] - rotstruct[i][1]*tempdiff[0];
      
        vecprod[0] = rotstruct[i][1]*xref[i][2] - rotstruct[i][2]*xref[i][1];
        vecprod[1] = rotstruct[i][2]*xref[i][0] - rotstruct[i][0]*xref[i][2];
        vecprod[2] = rotstruct[i][0]*xref[i][1] - rotstruct[i][1]*xref[i][0];
        
        
        for (j=0; j<3; j++) 
        {
            angmom[i][j] = mass[i]*vecprod[j];
        }
        
    }
    
    // make sum equal to zero
    for (i=0; i<3; i++) 
    {
        sum[i] = 0.0;
    }
    
    for (i=0; i<natom; i++) 
    {
        for (j=0; j<3; j++) 
        {
            sum[j] = sum[j] + angmom[i][j];
        }
    }
    
    for (k=0; k<3; k++) 
    {
        result = result + sum[k]*sum[k];
    }
    
    result = sqrt(result);
    */
        for (i=0; i<3; i++) 
        {
            sum[i] = 0.0;
        }
        
        for (i=0; i<natom; i++) 
        {
            vecprod[0] = -rotstruct[i][1]*xref[i][2] + rotstruct[i][2]*xref[i][1];
            vecprod[1] = -rotstruct[i][2]*xref[i][0] + rotstruct[i][0]*xref[i][2];
            vecprod[2] = -rotstruct[i][0]*xref[i][1] + rotstruct[i][1]*xref[i][0];
            
            
            for (j=0; j<3; j++) 
            {
                sum[j] = sum[j] + mass[i]*vecprod[j];
            }
            
        }
        
        for (i=0; i<3; i++) 
        {
            sum[i] = sum[i]*sum[i];
        }
        
        result = sqrt(sum[0] + sum[1] + sum[2]);
        
    return result;
    }
    
    //  numerical derivatives :
    void df(VecDoub &pokus, VecDoub &df)
    {
    /*
        // forward difference
        VecDoub x=pokus;
        Doub h=1e-4;
        Doub fold = operator()(pokus);
        
        for (Int i=0; i<4; i++) 
        {
            Doub temp = pokus[i];
            x[i] = temp + h;
            h = x[i] - temp;
            Doub fh = operator()(x);
            x[i] = temp;
            df[i] = (fh - fold)/h;
        }
       
      */  
        
       // central difference
        VecDoub x1 = pokus;
        VecDoub x2 = pokus;
        Doub h=1e-8;
        
        for (Int i=0; i<4; i++) 
        {
            Doub temp = pokus[i];
            x1[i] = temp + h;
            x2[i] = temp - h;
            h = x1[i] - temp;
            Doub fh1 = operator()(x1);
            Doub fh2 = operator()(x2);
            x1[i] = temp;
            x2[i] = temp;
            df[i] = (fh1 - fh2)/(2*h);
        }
        
        
        
    }
};
