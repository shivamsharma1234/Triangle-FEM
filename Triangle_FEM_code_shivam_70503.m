function []= Triangle_FEM_code_shivam_70503()
clc
%Shivam Sharma 70503
%Btech IT ( mathematical innovation)
start=input('Enter Start point of x & y');
n = input('End point of x: eg=3 ');
len2 = input('End point of y : ');
X = input('Triangular Elements Along X: ');
Y = input('Triangular Elements Along Y: ');
f0 = input('(F)not value:');


he  =  (n-start)/X;
ke = (len2-start)/Y;
node = (X+1)*(Y+1);
syms L1 L2 L3;
nodes = 2*(X*Y);



p = [0,0;  he,0; he,ke];
x1 = p(1, 1)*L1 + p(2, 1)*L2 + p(3, 1)*(1-L1-L2);
y1 = p(1, 2)*L1 + p(2, 2)*L2 + p(3, 2)*(1-L1-L2);

J = [diff(x1,L1), diff(y1,L1) ; diff(x1,L2) diff(y1,L2)];
Jinv=inv(J);
xi_diff_xy = inv(J)*[diff(L1, L1) diff(L2, L1) diff(1-L1-L2, L1); diff(L1, L2) diff(L2, L2) diff(1-L1-L2, L2)];

xi=[L1;L2;L3];
f0=0;
xi=[L1;L2;L3];

degree=input('Enter the degree')
for i=1:3
 	for j=1:3
        k_int(L1, L2, L3) = det(J)*(xi_diff_xy(1, i)*xi_diff_xy(1, j) + xi_diff_xy(2, i)*xi_diff_xy(2, j));
		K(i, j) = triangle_integrate(k_int,degree);
        
    end
    f_int(L1, L2, L3) = f0*xi(i)*det(J);
	F(i) = triangle_integrate(f_int);
end

disp('Matrix F single element:')
F
disp('Matrix K for a single element:')
K
F = F';
%conn=input('Enter the connectivity matrix')
%conn = [1 2 6; 2 3 7; 3 4 8; 1 6 5;2 7 6; 3 8 7; 5 6 10; 6 7 11;7 8 12;5 9 10;6 10 11;7 11 12];
conn = [1 2 6;1 5 6;2 3 7;2 6 7;3 4 8;3 7 8; 5 6 10;5 9 10;6 7 11;6 10 11;7 8 12;7 11 12];

disp('Assembled matrix for K')
K_assembled = Assemble_shivam_70503(K,conn)
disp('Assembled matrix for F')
F_assembled = Assemble_shivam_70503(F,conn)


u_index = input('Enter the boundary indexes of u whose value given : ');

b_q = input('Enter the Boundary value Q : ');

indi = u_index;
len = size(u_index);
len = len(2);
S2 = K_assembled;

for i=1:length(u_index)
    F_assembled=F_assembled-(K_assembled(:,u_index(i))*b_q(i));
end

F2 = F_assembled;



for i=1:len
    S2(u_index(i),:)=[];
    S2(:,u_index(i))=[];
    F2(u_index(i),:)=[];
    u_index = u_index-1;
end

K_decomposed = S2;

F_decomposed = F2;
U = (inv(K_decomposed)*(F_decomposed))

end



