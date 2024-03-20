clear all
close all
clc

A_a = 7e-4;
Q_a = 70.5;
C_a = 1e-38;
A_b = 10.1e-4;
Q_b = 14;
C_b = 4.26e-10;
R = 1.98720425864083e-3;
ReeErrying = @(dot_e,theta) A_a * (log(2*C_a*dot_e) + Q_a/(R*theta))  + ...
    A_b * asinh(C_b*dot_e*exp(Q_b/(R*theta)));


for numtemps = 3:6
    for numsamples = 3:6
        if numsamples*numtemps>18
            continue
        end

        temps = round(linspace(203,343,numtemps));

        xdata = logspace(-5,-2, numsamples);
        xdatafull = [];
        ydatafull = [];
        zdatafull = [];
        for j = 1:numel(temps)
            xdatafull = [xdatafull xdata];
            ydatafull = [ydatafull repmat(temps(j),size(xdata))];
            % zdatafull = [zdatafull ReeErrying(xdata,(temps(j))) + ReeErrying(xdata,(temps(j))).*(rand(size(xdata))*2-1)*0.1];
        end
        xdatafull = xdatafull';
        ydatafull = ydatafull';
        zdatafull = zdatafull';

        fprintf('Virtual experiment with %d samples and %d tempertures...',numsamples,numtemps);
        parfor i = 1:50000
            zdatafull = [];
            for j = 1:numel(temps)
                virtualError = ReeErrying(xdata,(temps(j))).*mean((rand(1,size(xdata,2))*2-1)*0.05,1);
                zdatafull = [zdatafull ReeErrying(xdata,(temps(j))) + virtualError];
            end
            zdatafull = zdatafull';


            [xData, yData, zData] = prepareSurfaceData( xdatafull, ydatafull, zdatafull );

            % Set up fittype and options.
            ft = fittype( 'A_a*7e-4 * (log(2*C_a*1e-38*x) + Q_a*70.5/(0.001987204258641*y)) + A_b*10.1e-4 * asinh(C_b*4.26e-10*x*exp(Q_b*14/(0.001987204258641*y)));', 'independent', {'x', 'y'}, 'dependent', 'z' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.DiffMinChange = 1e-16;
            opts.Display = 'Off';
            opts.Lower = [0 0 0 0 0 0];
            opts.Upper = ones(1,6)*1000;
            opts.MaxFunEvals = 600;
            opts.MaxIter = 400;
            opts.StartPoint = [1 1 1 1 1 1];
            opts.TolFun = 1e-12;
            opts.TolX = 1e-6;

            % Fit model to data.
            [fitresult, gof] = fit( [xData, yData], zData, ft, opts );

            Q_b = fitresult.Q_b*14;
            C_b = fitresult.C_b*4.26e-10;
            epsilon_crit = (exp(-Q_b/(R*273))/(2*C_b));
            error = (log10(epsilon_crit) - log10(0.007279483758))/log10(0.007279483758)*100;

            e_c_hist(i) = epsilon_crit;
            error_hist(i) = error;
        end

        fprintf('done\n')
        save(sprintf('error_hist_%d_sample_%d_temps.mat',numsamples,numel(temps)),'error_hist',"e_c_hist")
    end
end

%%
close all
clear all
clc
figure()
count = 1;
idx = [5 5 ; 3 4 ];
legendList={};
for i = 1:size(idx,1)

    legendList{count} =  sprintf('%d Tempertures and %d sample/temperature',idx(i,2),idx(i,1));
    fprintf('%d Tempertures and %d sample/temperature\n',idx(i,1),idx(i,2))
    count = count + 1;
    load(sprintf('error_hist_%d_sample_%d_temps.mat',idx(i,1),idx(i,2)))
    histogram(error_hist,'Normalization','probability')
    hold on

end
enhance_plot()
xlabel('Relative error  $\frac{log_{10}(\hat{\epsilon})-log_{10}(\epsilon)}{log_{10}(\epsilon)}$')
ylabel('Probability')
legend(legendList,'Location','north')
ylim([0 0.1])
title('Error Distribution')

%%
close all
figure()
count = 1;
idx = [5 5; 3 4 ];
legendList={};
Xmax = 0.2;
for i = 1:size(idx,1)
    
    legendList{count} =  sprintf('%d Tempertures and %d sample/temperature',idx(i,2),idx(i,1));
    count = count + 1;
    load(sprintf('error_hist_%d_sample_%d_temps.mat',idx(i,1),idx(i,2)))
    plotvar = log10(e_c_hist)-log10(0.007279483758);
    % plotvar = plotvar(plotvar<Xmax);
    histogram(plotvar,'Normalization','probability')
    hold on

end
enhance_plot()
xlabel('Deviation:  $log_{10}(\hat{\epsilon})-log_{10}(\epsilon)$')
ylabel('Probability')
legend(legendList,'Location','north')
ylim([0 0.09])
title('Error Distribution')
grid on
% %%
% close all
% figure()
% for i = 3:5
%     load(sprintf('error_hist_3_sample_%d_temps.mat',i))
%     histogram(error_hist,'Normalization','probability')
%     hold on
% end
% enhance_plot()
% xlabel('Relative error  $\frac{log_{10}(\hat{\epsilon})-log_{10}(\epsilon)}{log_{10}(\epsilon)}$')
% ylabel('Probability')
% legend('3 temperature','4 temperature', '5 temperature')
% title('3 samples')
% % temps
% % 3 temps Tempertures: 203 263 343
% % 4 temps Tempertures : 203 238 300 343
% %%
% close all
% figure()
% for i = 3:5
%     load(sprintf('error_hist_4_sample_%d_temps.mat',i))
%     histogram(error_hist,'Normalization','probability')
%     hold on
% end
% enhance_plot()
% xlabel('Relative error  $\frac{log_{10}(\hat{\epsilon})-log_{10}(\epsilon)}{log_{10}(\epsilon)}$')
% ylabel('Probability')
% legend('3 temperature','4 temperature', '5 temperature')
% title('4 samples')
% %%
% close all
% figure()
% for i = 3:5
%     load(sprintf('error_hist_5_sample_%d_temps.mat',i))
%     histogram(error_hist,'Normalization','probability')
%     hold on
% end
% enhance_plot()
% xlabel('Relative error  $\frac{log_{10}(\hat{\epsilon})-log_{10}(\epsilon)}{log_{10}(\epsilon)}$')
% ylabel('Probability')
% legend('3 temperature','4 temperature', '5 temperature')
% title('5 samples')
% %%
%
%%
close all

count = 1;
stdList = [];
for i = 3:6
    for j = 3:6
        try
            load(sprintf('error_hist_%d_sample_%d_temps.mat',i,j))
            stdList(i,j) = std(log10(e_c_hist));
        catch
            continue
        end
    end
end
disp(stdList(3:6,3:6))
%%




% % close all
% % y = 273;
% % figure()
% % x = logspace(-6,1);
% % semilogx(x,(fitresult.A_a*7e-4 * (log(2*fitresult.C_a*1e-38*x) + fitresult.Q_a*70.5/(0.001987204258641*y)) + fitresult.A_b*10.1e-4 * asinh(fitresult.C_b*4.26e-10*x*exp(fitresult.Q_b*14/(0.001987204258641*y))))*1e2,'b')
% % hold on
% % semilogx(x,(7e-4 * (log(2*1e-38*x) + 70.5/(0.001987204258641*y)) + 10.1e-4 * asinh(4.26e-10*x*exp(14/(0.001987204258641*y))))*1e2,color='black')
% % xline(0.007279483758)
% % xline(e_c_hist(end),Color='b')
%
%
% %%%
% % close all
% % figure()
% % for i = 1:numel(temps)
% %     semilogx(xdata,zdatafull(((i-1)*numsamples+1):i*numsamples)*1e2','.',Color='r')
% %     hold on
% %     x = xdata;
% %     y = temps(i);
% %     semilogx(x,(fitresult.A_a*7e-4 * (log(2*fitresult.C_a*1e-38*x) + fitresult.Q_a*70.5/(0.001987204258641*y)) + fitresult.A_b*10.1e-4 * asinh(fitresult.C_b*4.26e-10*x*exp(fitresult.Q_b*14/(0.001987204258641*y))))*1e2,'b')
% %     semilogx(x,(7e-4 * (log(2*1e-38*x) + 70.5/(0.001987204258641*y)) + 10.1e-4 * asinh(4.26e-10*x*exp(14/(0.001987204258641*y))))*1e2,color='black')
% % end
%
%
%
% % % Plot fit with data.
% % figure( 'Name', 'untitled fit 1' );
% % h = plot( fitresult, [xData, yData], zData );
% % legend( h, 'untitled fit 1', 'zdatafull vs. xdatafull, ydatafull', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % % Label axes
% % xlabel( 'xdatafull', 'Interpreter', 'none' );
% % ylabel( 'ydatafull', 'Interpreter', 'none' );
% % zlabel( 'zdatafull', 'Interpreter', 'none' );
% % grid on
% % view( -25.4, 52.1 );
%
%
%
%
%
%
%
%
%
%
%

