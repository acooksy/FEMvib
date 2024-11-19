//
//  no_reference.h
//  GMAT
//
//  Created by PETER ZAJAC on 4/2/13.
//  Copyright 2013 __MyCompanyName__. All rights reserved.
//



void trans_rot(MatDoub &xmat, MatDoub &interns, VecDoub &mass, Int &npt, Int &dim, Int &natom)

{
    
    // xmat     npt x 3*dim
    // mass     natom
    // interns  npt x dim
    
    
    Int i,j,k,l,m,h,pos;
    Int p,p_temp; // this will be used for indexing
    Int numstr=1; // this will always be 1 since we will go from point to point
    VecInt used(npt);
    VecInt order(npt);
    VecInt order_temp(npt);
    MatDoub ref(1,3*natom);
    MatDoub xgeom(1,3*natom);
    for (int q=0;q<npt;q++){used[q]=0;}
	Doub maxX, maxY, maxZ;
	Doub minX, minY, minZ;
    Doub radius, distance1, distance2, rr;
    Doub xcord, ycord, zcord;
    Int  how_many_pts=1;
    const double PI  = 3.141592653589793238462; 
    
    
    // case 3D:
    if (dim == 3)
    {
        MatDoub ns_i(1,3);
        maxX=interns[0][0];
        minX=interns[0][0];
        maxY=interns[0][1];
        minY=interns[0][1];
        maxZ=interns[0][2];
        minZ=interns[0][2];
        
        for (i=1; i<npt; i++) 
        {
            //check x
            if (interns[i][0] > maxX ){maxX = interns[i][0];}
            if (interns[i][0] < minX ){minX = interns[i][0];}
            
            //check y
            if (interns[i][1] > maxY ){maxY = interns[i][1];}
            if (interns[i][1] < minY ){minY = interns[i][1];}			
            
            //check z
            if (interns[i][2] > maxZ ){maxZ = interns[i][2];}
            if (interns[i][2] < minZ ){minZ = interns[i][2];}			
            
            // for loop end
        }   
        
        xcord = (maxX - minX);
        ycord = (maxY - minY);
        zcord = (maxZ - minZ);
        radius  =  pow(0.75*xcord*ycord*zcord/PI,1/3)/100.00; 
        
        VecDoub center_point(3);
        MatDoub start_pt(1,3);
        center_point[0] = minX + xcord/2;
        center_point[1] = minY + ycord/2;
        center_point[2] = minZ + zcord/2;
        
        howmany_noref(interns, center_point, dim, npt, how_many_pts, start_pt);
        // find order
        ordering(interns, order, start_pt, dim, npt, npt);
        
        used[order[0]]=1;
        for (i=1;i<npt; i++)  // not i=0 because that is the starting reference point
        {
            p=order[i];
            
            if (i==1)
            {
                // do the rotation for the first structure
                for (j=0; j<3*natom; j++) // define the reference structure
                {
                    ref[0][j] = xmat[order[0]][j];
                    xgeom[0][j] = xmat[p][j];
                }
                prerot(xgeom, ref, natom, numstr);
                rotation(ref, numstr , natom, xgeom, mass);
                playdata(xgeom, natom, numstr ,ref ); // ns_i - not important variable
                
                for (j=0; j<3*natom; j++) // define the reference structure
                {
                    xmat[p][j] = xgeom[0][j];
                }                
                used[p] = 1;
                
            }
            else
            {
                // search all the points again with the new one defined as reference:
                for (k=0; k<dim; k++) 
                {
                    start_pt[0][k] = interns[p][k];
                }
                
                ordering(interns, order_temp, start_pt, dim, npt, npt);
                
                // find which is the closest point that has been rotated already
                for (l=0; l<npt; l++) 
                {
                    p_temp = order_temp[l];
                    
                    if (used[p_temp] == 1)
                    {
                        for (j=0; j<3*natom; j++) // define the reference structure
                        {
                            ref[0][j] = xmat[p_temp][j];
                            xgeom[0][j] = xmat[p][j];
                        }
                        prerot(xgeom, ref, natom, numstr);
                        rotation(ref, numstr , natom, xgeom, mass);
                        playdata(xgeom, natom, numstr ,ref); // ns_i - not important variable
                        
                        for (j=0; j<3*natom; j++) // define the reference structure
                        {
                            xmat[p][j] = xgeom[0][j];
                        }                
                        used[p] = 1;
                        
                        break;
                        
                    }
                    
                    
                }
                cout << p_temp << " " << p << "\n";
                
            }
            
            
            
        }
        

        
        
            
    } // end 3D case
    
    
    // case 2D 
    if (dim == 2)
    {
        MatDoub ns_i(1,2);
        maxX=interns[0][0];
        minX=interns[0][0];
        maxY=interns[0][1];
        minY=interns[0][1];
        
        for (i=1; i<npt; i++) 
        {
            //check x
            if (interns[i][0] > maxX ){maxX = interns[i][0];}
            if (interns[i][0] < minX ){minX = interns[i][0];}
            
            //check y
            if (interns[i][1] > maxY ){maxY = interns[i][1];}
            if (interns[i][1] < minY ){minY = interns[i][1];}			
            
            // for loop end
        }
        
        xcord = (maxX - minX);
        ycord = (maxY - minY);
        radius  =  pow(xcord*ycord/PI,0.5)/100.00;
        
        VecDoub center_point(2);
        MatDoub start_pt(1,2);
        center_point[0] = minX + xcord/2;
        center_point[1] = minY + ycord/2;
        
        howmany_noref(interns, center_point, dim, npt, how_many_pts, start_pt);
        
     //   cout << "center point " << start_pt[0][0] << " " << start_pt[0][1] << "\n";
        
        
        // find order
        ordering(interns, order, start_pt, dim, npt, npt);
        
        used[order[0]]=1;
        for (i=1;i<npt; i++)  // not i=0 because that is the starting reference point
        {
            p=order[i];
            
            if (i==1)
            {
                // do the rotation for the first structure
                for (j=0; j<3*natom; j++) // define the reference structure
                {
                    ref[0][j] = xmat[order[0]][j];
                    xgeom[0][j] = xmat[p][j];
                }
                prerot(xgeom, ref, natom, numstr);
                rotation(ref, numstr , natom, xgeom, mass);
                playdata(xgeom, natom, numstr ,ref ); // ns_i - not important variable
             
                for (j=0; j<3*natom; j++) // define the reference structure
                {
                    xmat[p][j] = xgeom[0][j];
                }                
                used[p] = 1;
                
            }
            else
            {
                // search all the points again with the new one defined as reference:
                for (k=0; k<dim; k++) 
                {
                    start_pt[0][k] = interns[p][k];
                }
                
                ordering(interns, order_temp, start_pt, dim, npt, npt);
                
                // find which is the closest point that has already been rotated already
                for (l=0; l<npt; l++) 
                {
                    p_temp = order_temp[l];
                    
                    if (used[p_temp] == 1)
                    {
                        for (j=0; j<3*natom; j++) // define the reference structure
                        {
                            ref[0][j] = xmat[p_temp][j];
                            xgeom[0][j] = xmat[p][j];
                        }
                        prerot(xgeom, ref, natom, numstr);
                        rotation(ref, numstr , natom, xgeom, mass);
                        playdata(xgeom, natom, numstr ,ref ); // ns_i - not important variable
                        
                        for (j=0; j<3*natom; j++) // define the reference structure
                        {
                            xmat[p][j] = xgeom[0][j];
                        }                
                        used[p] = 1;
                        
                        break;
                        
                    }
                 
                    
                }
                cout << p_temp << " " << p << "\n";
                
            }
    
            
            
        }
        
                
    } // end case 2D    

    
    
    
    // end function 
}