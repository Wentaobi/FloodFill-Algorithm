clear,clc
map=zeros(1,6);
map(2:5,4)=-1;
y=[3,2];
map(3,2)=0;
z=[3,5];
 path=floodfillmex(map,y,z)