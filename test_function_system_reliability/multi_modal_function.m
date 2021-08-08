function obj=multi_modal_function(x)
obj1=-(x(:,1).^2+4).*(x(:,2)-1)/20+sin(5*x(:,1)/2)+2;
obj2=-(x(:,1)+2).^4+x(:,2)-4;
obj3=-(x(:,1)-4).^3+x(:,2)-2;

obj=[obj1,obj2,obj3];
end