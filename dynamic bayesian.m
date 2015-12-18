%model transition matrix

%list constant
state_variable=8;
node=4;
B=[1,-1,0,0,0,1,0,0;
   0,1,-1,0,0,0,1,-1;
   0,0,1,-1,0,-1,0,0;
   0,0,0,1,-1,0,-1,0]
%
a11=eye(state_variable);
a12=zeros(state_variable,node);
a21=B;
a22=eye(node)

A=[a11,a12;a21,a22]
print A