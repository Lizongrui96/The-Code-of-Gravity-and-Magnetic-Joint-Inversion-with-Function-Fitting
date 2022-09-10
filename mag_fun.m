function T=mag_fun(Vcc,xyz0,xyz,size0,IId,DDd,II0d,DD0d)
G=6.67*10^-2;
II = pi * IId / 180.;
DD = pi * DDd / 180.;
II0 = pi * II0d / 180.;
DD0 = pi * DD0d / 180.;
alphaM = cos(II) * cos(DD);
betaM = cos(II) * sin(DD);
gammaM = sin(II);
alphaT = cos(II0) * cos(DD0);
betaT = cos(II0) * sin(DD0);
gammaT = sin(II0);
%%
switch Vcc
    case 'Vt'
        vxx = grav_fun('Vxx',xyz0,xyz,size0); vxy = grav_fun('Vxy',xyz0,xyz,size0); 
        vxz = grav_fun('Vxz',xyz0,xyz,size0); vyy = grav_fun('Vyy',xyz0,xyz,size0); 
        vyz = grav_fun('Vyz',xyz0,xyz,size0); vzz = grav_fun('Vzz',xyz0,xyz,size0);
    case 'Vtx'
        vxx = grav_fun('Vxxx',xyz0,xyz,size0); vxy = grav_fun('Vxyx',xyz0,xyz,size0); 
        vxz = grav_fun('Vxzx',xyz0,xyz,size0); vyy = grav_fun('Vxyy',xyz0,xyz,size0); 
        vyz = grav_fun('Vxyz',xyz0,xyz,size0); vzz = grav_fun('Vxzz',xyz0,xyz,size0);
    case 'Vty'
        vxx = grav_fun('Vxyx',xyz0,xyz,size0); vxy = grav_fun('Vxyy',xyz0,xyz,size0); 
        vxz = grav_fun('Vxyz',xyz0,xyz,size0); vyy = grav_fun('Vyyy',xyz0,xyz,size0); 
        vyz = grav_fun('Vyzy',xyz0,xyz,size0); vzz = grav_fun('Vyzz',xyz0,xyz,size0);
    case 'Vtz'
        vxx = -grav_fun('Vxzx',xyz0,xyz,size0); vxy = -grav_fun('Vxyz',xyz0,xyz,size0); 
        vxz = -grav_fun('Vxzz',xyz0,xyz,size0); vyy = -grav_fun('Vyzy',xyz0,xyz,size0); 
        vyz = -grav_fun('Vyzz',xyz0,xyz,size0); vzz = -grav_fun('Vzzz',xyz0,xyz,size0);
end
bx = (alphaM*vxx + betaM*vxy + gammaM*vxz);
by = (alphaM*vxy + betaM*vyy + gammaM*vyz);
bz = (alphaM*vxz + betaM*vyz + gammaM*vzz);
dt = bx * alphaT + by * betaT + bz * gammaT;
T= dt * 100 / G;
end