function txyz=cross(Px,Py,Pz,CC,Px1,Py1,Pz1)
q=ones(size(CC'));
pyxc=(Py-Px)*CC*q;
pzyc=(Pz-Py)*CC*q;
pxzc=(Px-Pz)*CC*q;
txyz=pyxc.*Pz1+pzyc.*Px1+pxzc.*Py1;
end



