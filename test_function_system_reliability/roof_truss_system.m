function obj=roof_truss_system(x)   
         q=x(:,1);
         l=x(:,2);
         As=x(:,3);
         Ac=x(:,4);
         Es=x(:,5);
         Ec=x(:,6);
         fs=x(:,7);
         fc=x(:,8);
         obj(:,1)=0.03-0.5*q.*l.^2.*(3.81./(Ac.*Ec)+1.13./(As.*Es));
         obj(:,2)=1-1.185*q.*l./(fc.*Ac);
         obj(:,3)=1-0.75*q.*l./(fs.*As);

    
end