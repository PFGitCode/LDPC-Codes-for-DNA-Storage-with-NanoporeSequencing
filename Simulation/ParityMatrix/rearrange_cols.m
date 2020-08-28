function [b,rearranged_cols]=rearrange_cols(A)
%[b,rearranged_cols]=rearrange_cols(A)
%For examples and more details, please refer to the LDPC toolkit tutorial at
%http://arun-10.tripod.com/ldpc/ldpc.htm
%Rearrange the columns of the parity check matrix to get non singular A matrix

dim=size(A);
rows=dim(1);
cols=dim(2);
newA=A;
for i=1:rows
   rearranged_cols(i)=0;
end


   for i=1:rows
      if newA(i,i)==0
         k=i+1;
         if k<cols
         	%while A(rows,k)==0
             while newA(i,k)==0
               k=k+1;
            	if k==cols
               	break
            	end
         	end 
         	%if k~=cols
         		temp=newA(1:rows,k);
         		newA(1:rows,k)=newA(1:rows,i);
               newA(1:rows,i)=temp;
               
               rearranged_cols(i)=k;
               
               temp=A(1:rows,k);
         		A(1:rows,k)=A(1:rows,i);
            	A(1:rows,i)=temp;
         	%end
         end
      end
      for j=1:rows
        if j~=i
            if newA(j,i)==1
               newA(j,1:cols)=xor(newA(i,1:cols),newA(j,1:cols));
            end 
         end
     end
      A;

end   

rearranged_cols;
   b=A;