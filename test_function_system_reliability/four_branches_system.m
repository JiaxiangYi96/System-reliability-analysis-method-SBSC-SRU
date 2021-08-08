function obj=four_branches_system(x)
x1=x(:,1);
x2=x(:,2);
g1=3+(x1-x2).^2/10-(x1+x2)./sqrt(2);
g2=3+(x1-x2).^2/10+(x1+x2)./sqrt(2);
g3=(x1-x2)+6/sqrt(2);
g4=-(x1-x2)+6/sqrt(2);

g=[g1,g2,g3,g4];

obj1=min(g,[],2);
obj2=7*(x2+3).^2-5*x1.^2+(x1.^2+(x2+3).^2).^2+1;

obj3=2*(x1+2).^2-4*(x2-1)+((x1+2).^2+(x2-1).^2).^2+1;

obj=[obj1 obj2 obj3];

end