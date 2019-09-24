%D = dataset( 'file', './results_real_d/IQA/trueskill/20_5_15_16_19_26_28_31.csv', 'delimiter', ',' );
%D = dataset( 'file', './results_real_d/VQA/thurstone/1.csv', 'delimiter', ',' );
D = dataset( 'file', './results_rand_s/data/trueskill/200_0_5_16_19_28.csv', 'delimiter', ',' );
%pathSave = strcat('./results_rand_s/data/trueskill/nc_120_r_0_5_d_15_16_19_28');

n_exp = 12;

COLORs = lines( n_exp );
rangeysrocc = [0.5 0.8 0.9 0.95 0.98 0.995 0.999]; 
%rangeysrocc = [0.9 0.95 0.99 0.999 ];  %IQA
rangeyrmse =  [0.025 0.05 0.1 0.2 0.4 0.8]; 

LABELs = {};
ticks_style= {'-','--','-.',':'};
marker_style= {'>','o','*','s','+','x'};
ptypes = {'corr', 'rmse', 'mean_var' };
plabels = { 'SROCC', 'RMSE', 'Mean Variance' }; % [ -1/min(rangeyrmse) -1/max(rangeyrmse)]
ylimits = {[atanh(min(rangeysrocc)) atanh(max(rangeysrocc))],[log2(min(rangeyrmse)), log2(max(rangeyrmse))] };%_0_1

for pl=[1:2]
    H = figure(pl+2);
    html_change_figure_print_size( gcf, 16, 12 );
    clf
    pp = 1;
    for dd=1:n_exp
        Ds = D(D.design==dd,:);           
        if size(Ds,1)>0
            LT1 = ticks_style{mod(dd,size(ticks_style,2))+1};
            LT2 = marker_style{mod(dd,size(marker_style,2))+1};
            if( dd == 1 )
                LAB1 = 'FPC';
                marker_color = mean( [COLORs(dd,:); 1 1 1] );
            elseif dd == 2
                LAB1 = 'NC';
                marker_color = mean( [COLORs(dd,:); 1 1 1] );
            elseif dd == 3
                LAB1 = 'Swiss system';
                marker_color = 'none';
            elseif dd == 4
                LAB1 = 'Adaptive squares';
                marker_color = 'none';
            elseif dd == 5
                LAB1 = 'Peng Ye';
                marker_color = 'none';
            elseif dd == 6
                LAB1 = 'Quicksort';
                marker_color = 'none';
            elseif (dd==7 || dd==31)
                LAB1 = 'TS-sampling';
                marker_color = 'none';
            elseif (dd==8 ||dd==15)
                LAB1 = 'Crowd-BT';
                marker_color = 'none';
            elseif (dd==9 ||dd==16)
                LAB1 = 'HR-active';
                marker_color = 'none';
            elseif (dd==10 || dd == 19)
                LAB1 = 'Hybrid-MST';
                marker_color = 'none';  
            elseif (dd==11 || dd == 26)
                LAB1 = 'ASAP';
                LT2 = 'd';
                marker_color =mean( [COLORs(dd,:); 1 1 1] );
            elseif (dd==12 || dd==28)
                LT2 = 'd';
                LAB1 = 'ASAP-approx';% KL online true';
                marker_color = 'none';
            elseif (dd==13)
                LT2 = 'd';
                LAB1 = 'ASAP without selective EIG';% KL online true';
                marker_color = 'none';
            end
            
            if pl==1
                yyplot = atanh(Ds.(ptypes{pl}));
            elseif pl == 2
                yyplot = log2(Ds.(ptypes{pl}));% -1./(Ds.(ptypes{pl}));
            else
                yyplot =  -1./Ds.(ptypes{pl});
            end
            
            plot( Ds.cmps_per_n_conds, yyplot, 'LineStyle',LT1 ,'Marker', LT2, 'Color', COLORs(dd,:), 'MarkerFaceColor', marker_color   );
            hold on;
            LABELs{pp} =  LAB1;
            
            %for xx = 0:1:15
            %   xxlabel{xx+1} = num2str(xx);
            %end
            %xlim([0 15]);
            %set(gca, 'XTick', [0:1:15]), 
            %set(gca,'XTickLabel', xxlabel),
            

            
            %set(gca, 'YScale', 'log')
            pp = pp+1;
        end
        
    end
        d=0;
        yylabel= {};
        if pl == 1
            for yy = rangeysrocc
                d = d+1;
                yylabel{d} = num2str(yy);
            end

            set(gca,'YTick', atanh(rangeysrocc)),  
            set(gca,'YTickLabel',yylabel);

        elseif pl ==2
            for yy = rangeyrmse
                d = d+1;
                yylabel{d} = num2str(yy);
            end
            set(gca,'YTick', log2(rangeyrmse)),
            set(gca,'YTickLabel',yylabel); 
        else
           for yy = rangeyvar
                d = d+1;
                yylabel{d} = num2str(yy);
            end
            set(gca,'YTick', -1./rangeyvar),  
            set(gca,'YTickLabel',yylabel); 
            ;
        end
    grid on
    %xticks([0.5:0.25:2.5])
    %xticklabels({'0.5','0.75','1','1.25','1.5','1.75','2.0','2.25','2.5'})
    %xticks(1:15)
    %xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'})
    xticks(1:16)
    xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'})
    %xticks([0.1:0.1:1])
    %xticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
    xlabel( 'Standard trial');
    ylabel( plabels{pl});
    legend( LABELs, 'Location', 'best');
    ylim( ylimits{pl} );
    %xlim([0 15])
    %xlim([0.1 1])
    set(findall(H,'-property','FontSize'),'FontSize',16)
    %pathSave_n = strcat(pathSave,plabels{pl});
    %saveas(gcf,[pathSave_n],'jpg');
    
end