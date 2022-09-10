function Vij = grav_fun(Vcc, xyz0,xyz,size0)
G=6.67*10^-6;
o=1e-3;
%% 观测点与立方体最小点与最大点的相对位置
obs_model=[xyz0-xyz+size0./2;xyz0-xyz-size0./2]';
%%
Vij = 0;
for i = 1:2
    x = obs_model(1,i);
    x2 = x*x;
    for j = 1:2
        y = obs_model(2,j);
        y2 = y*y;
        for k = 1:2
            z = -(obs_model(3,k));
            z2 = z*z;
            r = sqrt(x2 + y2 +z2);
            uijk = (-1)^(i + j + k );
            switch Vcc
                case 'Vx'
                    if(x==0 && y==0 && z==0)
                        x=o;y=o;z=o;r=sqrt(x*x+y*y+z*z);
                    elseif(x==0)
                        x=o;r=sqrt(x*x+y*y+z*z);
                    end
                    if(y+r<=0)
                        z=o;r=sqrt(x*x+y*y+z*z);
                    elseif(z+r<=0)
                        y=o;r=sqrt(x*x+y*y+z*z);
                    end
                    Vtemp=z*log(y+r)+y*log(z+r)-x*atan(z*y/(x*r));
                case 'Vy'
                    if(x==0 && y==0 && z==0)
                        x=o;y=o;z=o;r=sqrt(x*x+y*y+z*z);
                    elseif(y==0)
                        y=o;r=sqrt(x*x+y*y+z*z);
                    end
                    if(x+r<=0)
                        z=o;r=sqrt(x*x+y*y+z*z);
                    elseif(z+r<=0)
                        x=o;r=sqrt(x*x+y*y+z*z);
                    end
                    Vtemp=z*log(x+r)+x*log(z+r)-y*atan(z*x/(y*r));
                case 'Vz'
                    if(x==0 && y==0 && z==0)
                        x=o;y=o;z=o;r=sqrt(x*x+y*y+z*z);
                    elseif(z==0)
                        z=o;r=sqrt(x*x+y*y+z*z);
                    end
                    if(x+r<=0)
                        y=o;r=sqrt(x*x+y*y+z*z);
                    elseif(y+r<=0)
                        x=o;r=sqrt(x*x+y*y+z*z);
                    end
                    Vtemp=y*log(x+r)+x*log(y+r)-z*atan(x*y/(z*r));
                case 'Vxx'
                    if(x==0)
                        Vtemp=-1e4*pi/2;
                    else
                        Vtemp=-1e4*atan(y*z/(x*r));
                    end
                case 'Vyy'
                    if(y==0)
                        Vtemp=-1e4*pi/2;
                    else
                        Vtemp=-1e4*atan(x*z/(y*r));
                    end
                case 'Vzz'
                    if (z==0)
                        Vtemp=-1e4*pi/2;
                    else
                        Vtemp=-1e4*atan(x*y/(z*r));
                    end
                case 'Vxz'
                    if (y+r==0)
                        x=o;r=sqrt(x*x+y*y+z*z);
                    end
                    Vtemp=1e4*log(y+r);
                case 'Vyz'
                    if (x+r==0)
                        y=o;r=sqrt(x*x+y*y+z*z);
                    end
                    Vtemp=1e4*log(x+r);
                case 'Vxy'
                    if (z+r==0)
                        x=o;r=sqrt(x*x+y*y+z*z);
                    end
                    Vtemp=1e4*log(z+r);
                case 'Vxxx'
                    if(x==0 && y==0)
                        x=o;y=o;
                    end
                    if(z==0 && x==0)
                        z=o;
                    end
                    r=sqrt(x*x+y*y+z*z);
                    Vtemp=z*y*(r*r+x*x)/((z*z+x*x)*(y*y+x*x)*r);
                case 'Vxyx'
                    if(y==0 && x==0)
                        y=o;
                    end
                    r=sqrt(x*x+y*y+z*z);
                    Vtemp=-z*x/((x*x+y*y)*r);
                case 'Vxzx'
                    if(x==0 && z==0)
                        z=o;
                    end
                    r=sqrt(x*x+y*y+z*z);
                    Vtemp=-x*y/((x*x+z*z)*r);
                case 'Vxyy'
                    if(y==0 && x==0)
                        x=o;
                    end
                    r=sqrt(x*x+y*y+z*z);
                    Vtemp=-z*y/((y*y+x*x)*r);
                case 'Vyyy'
                    if(x==0 && y==0)
                        x=o;y=o;
                    end
                    if(z==0 && y==0)
                        z=o;
                    end
                    r=sqrt(x*x+y*y+z*z);
                    Vtemp=x*z*(r*r+y*y)/((x*x+y*y)*(y*y+z*z)*r);
                case 'Vyzy'
                    if(y==0 && z==0)
                        z=o;
                    end
                    r=sqrt(x*x+y*y+z*z);
                    Vtemp=-x*y/((y*y+z*z)*r);
                case 'Vxzz'
                    if(x==0 && z==0)
                        z=o;
                    end
                    r=sqrt(x*x+y*y+z*z);
                    Vtemp=y*z/((z*z+x*x)*r);
                case 'Vyzz'
                    if(y==0 && z==0)
                        z=o;
                    end
                    r=sqrt(x*x+y*y+z*z);
                    Vtemp=x*z/((z*z+y*y)*r);
                case 'Vzzz'
                    if((x==0 || y==0) && z==0)
                        z=o;r=sqrt(x*x+y*y+z*z);
                    end
                    Vtemp=x*y*(r*r+z*z)/((x*x+z*z)*(y*y+z*z)*r);
                case 'Vxyz'
                    if(x==0 && y==0 && z==0)
                        r=o;
                    end
                    Vtemp=1/r;
            end
            Vij = Vij + uijk*Vtemp;
        end
    end
end
Vij = G*Vij;
end