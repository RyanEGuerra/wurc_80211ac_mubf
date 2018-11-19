function [] = mySaveAs(hand, savepath, width, height, varargin)

% PATH SHOULD NOT HAVE FILE EXTENSION!!

% figure(hand)
% 
% 
% set(gca,'units','inches')
% pos = get(gca,'Position');
% ti = get(gca,'TightInset');
% 
% set(gcf, 'PaperUnits','inches');
% set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)]);%+ti(2)+ti(4)]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)]);%+ti(2)+ti(4)]);

if(length(varargin)>0)
   EN_DEBUG = varargin{1};
    EN_DEBUG
else
    EN_DEBUG = 0;
end




set(hand, 'PaperUnits', 'inches');
set(hand, 'Units', 'inches');

ti = get(gca, 'TightInset');




set(hand, 'PaperSize', [width height]);

% if(EN_DEBUG==0)
set(hand, 'Position',[0 0 width height])
% end

set(hand, 'PaperPosition',[0 0 width height]);

if(EN_DEBUG==0)
% epspath = [savepath '.eps'];
pdfpath = [savepath '.pdf'];
saveas(hand,pdfpath, 'pdf') ;
end
%  saveas(hand,[savepath '.eps'], 'epsc2') ;
%  system(['epstopdf ' epspath]);
 
 
 
%  saveas(hand,[savepath '.pdf'], 'pdf') ;
 

% [status, result] = system(cmd);
 %disp(['i_view32 ' epspath ' /convert=' pdfpath])

%saveas(hand,[savepath '.emf'], 'emf') ;
 
end