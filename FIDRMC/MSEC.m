function result = MSEC(A,B,rand);
[mA nA,d] = size(A);
[mB nB,d] = size(B);
if (mA ~= mB | nA ~= nB)
   error('Arrays A and B must have the same size');
end;
if (d ~= 3)
   error('Only colour images!');
end;
A = double(A);
B = double(B);
result1 = 0;
result2 = 0;
result3 = 0;
% we laten de rand buiten beschouwing
for i = 1+rand:mA-rand
   for j = 1+rand:nA-rand
      result1 = result1 + (A(i,j,1) - B(i,j,1))^2;
      result2 = result2 + (A(i,j,2) - B(i,j,2))^2;
      result3 = result3 + (A(i,j,3) - B(i,j,3))^2;
   end;
end;
result1 = result1 / ((mA-2*rand)*(nA-2*rand));
result2 = result2 / ((mA-2*rand)*(nA-2*rand));
result3 = result3 / ((mA-2*rand)*(nA-2*rand));
result = (result1+result2+result3)/3.0;