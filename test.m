Active_GE_pos=zeros(1,4);                                   
%chose which instrrument has the active gross error         

%----------|
pos=[1 2 3 ]';   %|
%----------|
sum=0;
for i=1:4
    sum=Active_GE_pos(pos(j));
if (~isempty(pos))
    for j=1:length(pos)

       Active_GE_pos(pos(j))=1;
       
    end  
end
end
Active_GE_pos

sum