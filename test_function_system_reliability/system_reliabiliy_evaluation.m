function [pf_real,cov_estimate]=system_reliabiliy_evaluation(search_x,test_function,type)
% this function is used to evaluate the performance of the current
% iteration
 switch type
     case 'parallel'
    num_search=size(search_x,1);
    %     for ii=1:num_obj
    %         search_y(:,ii)=predictor(search_x,kriging_obj{ii});
    %     end
    %     num_fail=sum(sum(search_y < 0, 2) == size(search_y,2));
    real_y=feval(test_function,search_x);
    num_real = sum(sum(real_y < 0, 2) == size(real_y,2));
    %     pf_estimate=num_fail/num_search;
    pf_real=num_real/num_search;
        cov_estimate=sqrt((1-pf_real)/(num_search*pf_real));
    %     real_relative_error=abs(pf_estimate-pf_real)./pf_real;
     case  'series'
    num_search=size(search_x,1);
    %     for ii=1:num_obj
    %         search_y(:,ii)=predictor(search_x,kriging_obj{ii});
    %     end
    %     num_fail=sum(sum(search_y > 0, 2)> 0);
    real_y=feval(test_function,search_x);
    num_real = sum(sum(real_y < 0, 2)>0);
    %     pf_estimate=num_fail/num_search;
    pf_real=num_real/num_search;
        cov_estimate=sqrt((1-pf_real)/(num_search*pf_real));
    %     real_relative_error=abs(pf_estimate-pf_real)./pf_real;
 end
end