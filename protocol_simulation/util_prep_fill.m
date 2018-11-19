function [X, Y] = prep_fill(x, y1, y2)
%     
%     y1 = yu(:)+ys(:);
%     y2 = yu(:)-ys(:);
    
    X = [x(:); flipud(x(:))];
    Y = [y1(:); flipud(y2(:))];
    
    
    
    
    
end