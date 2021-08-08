function[pf_estimate,pf_real,cov_estimate,real_relative_error]=reliabiliy_evaluation_single_fidelity(search_x,kriging_model,test_function)
% this function is used to evaluate the performance of the current
% iteration
num_search=size(search_x,1);
search_y=predictor(search_x,kriging_model);
real_y=feval(test_function,search_x);
num_fail=sum(search_y<0);
num_real=sum(real_y<0);
pf_estimate=num_fail/num_search;
pf_real=num_real/num_search;
cov_estimate=sqrt((1-pf_estimate)/(num_search*pf_estimate));
real_relative_error=abs(pf_estimate-pf_real)./pf_real;
end