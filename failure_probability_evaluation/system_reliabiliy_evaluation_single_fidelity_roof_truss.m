function [pf_estimate,pf_real,cov_estimate,real_relative_error]=system_reliabiliy_evaluation_single_fidelity_roof_truss(search_x,num_obj,kriging_obj,test_function,type)
% this function is used to evaluate the performance of the current
% iteration
switch type
    case 'parallel'
        num_search=size(search_x,1);
        for ii=1:num_obj
            search_y(:,ii)=predictor(search_x,kriging_obj{ii});
        end
        num_fail=sum(sum(search_y < 0, 2) == size(search_y,2));
        real_y=feval(test_function,search_x);
        num_real = sum(sum(real_y < 0, 2) == size(real_y,2));
        pf_estimate=num_fail/num_search;
        pf_real=num_real/num_search;
        cov_estimate=sqrt((1-pf_estimate)/(num_search*pf_estimate));
        real_relative_error=abs(pf_estimate-pf_real)./pf_real;
    case 'series'
         num_search=size(search_x,1);
         search_x_1=search_x(:,1:6);
         search_x_2=[search_x(:,1:2),search_x(:,4),search_x(:,8)];
         search_x_3=[search_x(:,1:2),search_x(:,3),search_x(:,7)];          
         search_y(:,1)=predictor(search_x_1,kriging_obj{1});
         search_y(:,2)=predictor(search_x_2,kriging_obj{2});
         search_y(:,3)=predictor(search_x_3,kriging_obj{3});
          num_fail=sum(sum(search_y < 0, 2)> 0);
          real_y=feval(test_function,search_x);
          num_real = sum(sum(real_y < 0, 2) >0);
         pf_estimate=num_fail/num_search;
         pf_real=num_real/num_search;
         cov_estimate=sqrt((1-pf_estimate)/(num_search*pf_estimate));
         real_relative_error=abs(pf_estimate-pf_real)./pf_real;
        
end
end