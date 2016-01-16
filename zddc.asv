function [hl,t]=zddc(A)
count=1;
[n n] =size(A);

for p=1:n

h(p)=det(A(1:p, 1:p));

end

hl=h(1:n);zA=A';

for i=1:n

   if h(1,i)<=0 
    
    disp('pay attention: det(Ai) not all is >0, A is not'), hl;zA,return
    
end

end

if h(1,i)>0

disp('A is positive. transptation of A and leading positive is?')
t=count;
hl;zA
t;
end
end

