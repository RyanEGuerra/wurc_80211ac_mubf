function [ ] = plotthis( varargin )
%function [ ] = plot_these_vectors( varargin )
%   arg 1 = vector array of arbitrary size; I'll figure it out
%   arg 2 = optional label string to append to each plot title
%
% (c) ryan@guerra.rocks 2015
% http://www.apache.org/licenses/LICENSE-2.0

% Vectors to plot
X = varargin{1};

% addtional title for plotting
if nargin > 1
    suffix = varargin{2};
else
    suffix = '';
end

% debug
if nargin > 2
    warning([mfilename ': Only using first two arguments.']);
end

%
if isempty(X)
    error([mfilename ': Input is empty']);
end

% find the dimention of X that is the time sample dimension (to handle
% transpositions, etc...
if length(size(X)) > 2
    error('Input matrix is 3D... No can do.')
end
dim = find(size(X)==max(size(X)));
if dim==2
    X = transpose(X);
end


figure()
    numplots = size(X, 2);
    for ii = 1:1:numplots
        ax(ii) = subplot(numplots, 1, ii);
            if myisreal(X(:, ii))
                % kill any residual/rounding imaginary component
                plot(real(X(:, ii)));
            else 
                % we treat it as a complex vector
                plot(abs(X(:, ii)));
            end
            title(['Vector: ' num2str(ii) ' ' suffix]);
            grid on;
    end
    
    linkaxes(ax, 'x');
end

function [ its_real ] = myisreal( in_vecs )
%function [ its_real ] = myisreal( in_vecs )
%
%   Tests an input array to tell if there exist imaginary values within it,
%   to some tolerance.

    tol = 1e-4;
    res = abs(imag(in_vecs)) > tol;
    if sum(res(:)) > 0
        its_real = false;
    else
        its_real = true;
    end
end