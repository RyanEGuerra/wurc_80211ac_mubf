function [] = mySaveAs(hand, savepath, width, height, varargin)
% Use 3 and 5 for width and height with default matlab figure font sizes

if(isempty(varargin))
    debug = 0;
else
    debug = 1;
    disp(' -- SAVE AS DEBUG MODE --')
end

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


set(hand, 'PaperUnits', 'inches');
set(hand, 'Units', 'inches');

ti = get(gca, 'TightInset');




set(hand, 'PaperSize', [width height]);
set(hand, 'Position',[0 0 width height])
set(hand, 'PaperPosition',[0 0 width height]);


% epspath = [savepath '.eps'];
pdfpath = [savepath '.pdf'];

% if(debug==0)
%  saveas(hand,[savepath '.eps'], 'epsc2') ;
%  system(['epstopdf ' epspath]);
% end
 if(debug==0)
 
 saveas(hand,[savepath '.pdf'], 'pdf') ;
 end

% [status, result] = system(cmd);
 %disp(['i_view32 ' epspath ' /convert=' pdfpath])

%saveas(hand,[savepath '.emf'], 'emf') ;
 
end