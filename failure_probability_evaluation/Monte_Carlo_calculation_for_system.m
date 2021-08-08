clc
clear all
%% import basic information of the test function
test_function='series_three_modes';
[num_vari,mu,sigma,design_space,type]=test_function_system_reliability(test_function);
%% Calculate the failure probability
num_search=10^6;
search_x=MCS_Population_Generation(mu,sigma,num_search);
pf_real=system_reliabiliy_evaluation_single_fidelity(search_x,test_function,type);

%%

figure (1)
[x1_plot,x2_plot]=meshgrid(-5:0.005:5,-5:0.005:5);
ind=(8*x2_plot.^2-8*x1_plot.^2+(x1_plot.^2+x2_plot.^2).^2<0) & (2*x1_plot.^2-2*x2_plot.^2-(x1_plot.^2+x2_plot.^2).^2<0)...
    &(8*x2_plot.^2-8*x1_plot.^2-(x1_plot.^2+x2_plot.^2).^2<0);
plot(x1_plot(ind),x2_plot(ind),'.','MarkerEdgeColor',[0.3010 0.7450 0.9330],'MarkerSize',2);
hold on 

num_test=100;
x1_plot=linspace(design_space(1,1),design_space(2,1),num_test);
x2_plot=linspace(design_space(1,2),design_space(2,2),num_test);
% calculate the real responses
for ii=1:1:num_test
    for jj=1:1:num_test
        X_temp=[x1_plot(ii),x2_plot(jj)];
        Y_temp=feval(test_function,X_temp);
        y_real_1(jj,ii)=Y_temp(:,1);
        y_real_2(jj,ii)=Y_temp(:,2);
        y_real_3(jj,ii)=Y_temp(:,3);
    end
end



figure (1)
v=[-1e6,0,1e6];
contour(x1_plot,x2_plot,y_real_1,v,'k--','LineWidth',2,'ShowText','on');
hold on
contour(x1_plot,x2_plot,y_real_2,v,'k--','LineWidth',2,'ShowText','on');
hold on
contour(x1_plot,x2_plot,y_real_3,v,'k--','LineWidth',2,'ShowText','on');
hold on
set(gca,'XTick',design_space(1,1):(design_space(2,1)-design_space(1,1))/5:design_space(2,1))
set(gca,'YTick',design_space(1,2):(design_space(2,2)-design_space(1,2))/5:design_space(2,2))
set(gca,'fontname','Times New Roman','LineWidth',1.5,'fontsize',18)
xlabel('\it x_1','fontname','Times New Roman','fontsize',16)
ylabel('\it x_2','fontname','Times New Roman','fontsize',16)
