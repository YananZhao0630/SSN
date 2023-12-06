function [FA,MA,TE,Pfa,Pma,Perr,Pdet,DA,Kappa,OverallAccuracy,UserAccuracy,ProducerAccuracy] = counterror(I,M)

%  counterror
%  
%  [FA,MA,TE,Pfa,Pma,Perr,Pdet,DA,Kappa,OverallAccuracy,UserAccuracy,...    
%  ProducerAccuracy] = counterror(I,M, filename, path)
%  
%  Input
%      I - Binary image (logical Matlab data type)
%      M - Binary image taken as reference (logical Matlab data type)
%      filename - name of the file
%      path - path where the results will be saved
%  
%  Output
%      FA - Number of false alarms
%      MA - Number of missed alarms
%      TE - Number of total errors (FA+MA)
%      Pfa - Probability of false alarms
%      Pma - Probability of missed alarms
%      Perr - Probability of error
%      Pdet - Probability of detection
%      DA - Number of detected alarms
%      Kappa - Kappa accuracy
%      OverallAccuracy - Overall Accuracy
%      UserAccuracy - User Accuracy (2 values)
%      ProducerAccuracy - Producer Accuracy (2 values)
%  
%   Function that calculates the number of
%   false alarms (FA), missed alarms (MA),
%   total errors (TE), the probability of error (Perr),
%   the probability of detection (Pdet) and the detected alarms (DA)
%   comparing the change detection binary map I with the reference map M.
%
%           WARNING:
%           Change   0
%           Unchange 1
%  
%  -----------------
%  Mauro Dalla Mura
%  19 Jan 2007

[row col] = size(I);

CH = 0; % # of changed pixels (according to M, the reference map taken as ground truth)
FA = 0; % # of false alarms
MA = 0; % # of missed alarms
DA = 0; % # of detected alarms (pixels correctly classified to the changed class)
TE = 0; % # of total errors

% Mrgb(:,:,1) = uint8(I);
% Mrgb(:,:,2) = uint8(I);
% Mrgb(:,:,3) = uint8(I);

for i=1:row
    for j=1:col
        if(I(i,j)==1)&&(M(i,j)==0)
            FA = FA+1;
%             pix = [255 0 0];    % red
%             Mrgb(i,j,:)=pix;
        end
        if(I(i,j)==0)&&(M(i,j)==1)
            MA = MA+1;
            CH = CH+1;
%             pix = [0 255 0];    % green
%             Mrgb(i,j,:)=pix;        
        end
%         if(I(i,j)==1)&&(M(i,j)==1)
%             pix = [255 255 255];    % white
%             Mrgb(i,j,:)=pix;                      
%         end
        if(I(i,j)==1)&&(M(i,j)==1)
            DA = DA+1;
            CH = CH+1;
%             pix = [0 0 255];    % blue
%             Mrgb(i,j,:)=pix;                      
        end
    end
end

TE = FA+MA;
Perr = TE*100/(row*col);
Pdet = DA*100/CH;
Pfa = FA*100/(row*col-CH);
Pma = MA*100/CH;

% Compute Kappa accuracy
UNCH = row*col - CH;    % # of unchanged pixels (according to M, the reference map taken as ground truth)

m11 = UNCH - FA;
m12 = MA;
m21 = FA;
m22 = CH - MA;

r1 = m11 + m12;
r2 = m21 + m22;
c1 = m11 + m21;
c2 = m12 + m22;

UserAccuracy(1) = m11*100/r1;
UserAccuracy(2) = m22*100/r2;
ProducerAccuracy(1) = m11*100/c1;
ProducerAccuracy(2) = m22*100/c2;

a = m11 + m22;
b = c1*r1 + c2*r2;
a1 = a/(r1 + r2);
b1 = b/(r1 + r2)^2;

OverallAccuracy = a1*100;
Kappa = (a1 - b1)/(1 - b1);


% f = char(strcat([path,filename],'_ErrMap'));
% imwrite(Mrgb,[f,'.tif'],'tif');
% figure
% iptsetpref('ImshowBorder','tight');
% imshow(Mrgb);
% iptsetpref('ImshowBorder','loose');
% saveas(gcf,[f,'.eps'], 'psc2');
% close;
% 
% figure1 = figure(...
%   'PaperPosition',[0 0 20 20],...
%   'PaperSize',[20.98 29.68]);
% %'PaperPosition',[0 0 20 20],... for 800x800 img
% iptsetpref('ImshowBorder','tight');
% imshow(Mrgb);
% hold on
% bar(1,0,'b');
% bar_handle=bar(1,0,'r');
% bar(1,0,'g');
% baseline_handle = get(bar_handle,'BaseLine');
% set(baseline_handle,'LineStyle','none');
% legend1 = legend(['DA ',int2str(DA)],['FA ',int2str(FA)],['MA ',int2str(MA)],'Location','NorthWest');
% set(legend1, 'Box', 'off');
% 
% %% Create textbox
% annotation1 = annotation(...
%   figure1,'textbox',...
%   'Position',[0.05 0.832 0.115 0.0525],...
%   'FitHeightToText','off',...
%   'FontWeight','bold',...
%   'LineStyle','none',...
%   'String',{['Pdet ',num2str(Pdet,'%3.3f%%')],['Perr ',num2str(Perr,'%3.3f%%')]});
%  %  'Position',[0.03 0.8638 0.115 0.0525],...
% 
% %% Create rectangle
% %annotation2 = annotation(figure1,'rectangle',[0.0075 0.8638 0.17 0.1275]);
% annotation2 = annotation(figure1,'rectangle',[0.02 0.827 0.2 0.16]);
% 
% hold off
% iptsetpref('ImshowBorder','loose');
% saveas(gcf,[path,'I_ErrMap_',filename,'.eps'], 'psc2');
% saveas(gcf,[path,'I_ErrMap_',filename,'.tif'], 'tif');
% close;