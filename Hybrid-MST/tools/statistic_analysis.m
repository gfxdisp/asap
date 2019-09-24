function [CC, MAE, RMSE, ROCC]=statistic_analysis(predicted,original)                 
                      a=predicted;
                      b=original;
%                    a=imquality_i';
                   % b=100-dmos_i';
%                     b=dmos_i';
                    %%%%%%%%%%%%%%%
                    fun=inline('(beta(1)-beta(2))./(1+exp(-(a-beta(3))/abs(beta(4))))+beta(2)','beta','a');
                    beta0=[max(b) min(b) sum(a)/length(a) 1]';
                    [betap,R,J] = nlinfit(a,b,fun,beta0);
                    [DMOSp,delta] = nlpredci(fun,a,betap,R,J);
                    M_dmos_i = sum(b)/length(b);
                    M_DMOSp = sum(DMOSp)/length(DMOSp);
                    CC=sum((DMOSp'-M_DMOSp).*(b'-M_dmos_i))/sqrt(sum((DMOSp'-M_DMOSp).^2)*sum((b'-M_dmos_i).^2));
                    MAE=sum(abs(DMOSp-b))/size(b',2);
                    RMSE=sqrt(sum((DMOSp-b).^2)/size(b',2));
                    %OR(i)=size(find(res>2+std()b'),2)/size(b',2);
                    [imquality_i_ord oe_ord1] =sort(a);
                    [res oe_ord]=sort(oe_ord1);
                    [dmos_i_ord se_ord1]=sort(b);
                    [res se_ord]=sort(se_ord1);
                    ROCC=1-6*sum((oe_ord-se_ord).^2)/(length(oe_ord)*(length(oe_ord)^2-1)) ;
                    
%                     figure;
%                     plot(a,b,'k+');
%                     hold on
%                     plot(imquality_i_ord,DMOSp(oe_ord1),'k-');
%                     hold on
%                     xlabel('predicted'),ylabel('original')