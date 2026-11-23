function ModeshapeAnimation(x3, modeshape_disp, cp, xi, options)
%Make an animation for a modeshape
% Note: pure SH waves cannot be animated through this function because in pure SH waves, modeshape_u1=0 and modeshape_u3=0,
% which further causes the variable factor_scale computed in line 71 to be Inf.
% Author: Shuanglin GUO @PIMM of ENSAM
% May 8, 2020
% Email: shuanglin.guo@ensam.eu

% set default parameters for options---------------------------------------
options_default.num_wavelength = 2;      % num_wavelength must be a positive integer
options_default.num_period = 2;          % num_period must be a positive integer
options_default.scale = 1;               % scale must be a positive number
options_default.nondeformedgrids = 'ON'; % 2 options: 'ON' and 'OFF'
options_default.savegif = 'OFF';         % 2 options: 'ON' and 'OFF'
options_default.filename = 'untitled_ModeshapeAnimation2D.gif';  % filename to save the gif file
if nargin<5
    options = options_default;
else
    options = setOptions(options, options_default);
end
num_wavelength = options.num_wavelength;
num_period = options.num_period;
scale = options.scale;
nondeformedgrids = options.nondeformedgrids;
savegif = options.savegif;
filename = options.filename;
% end to set default parameters for options--------------------------------


if ~endsWith(filename,'.gif')
    filename = [filename, '.gif'];
end

% algorithm begins here----------------------------------------------------
fq = cp*xi/(2*pi);
wavelen = 2*pi/xi;  % wavelength with unit: m

if num_wavelength==1
    x1 = linspace(0, num_wavelength*wavelen, 25);
else
    x1 = linspace(0, num_wavelength*wavelen, 25*num_wavelength-1);
end
if num_period==1
    Time = linspace(0, num_period/fq, 13);
else
    Time = linspace(0, num_period/fq, 13*num_period-1);
end

num_x1 = length(x1); num_x3 = length(x3);

if size(modeshape_disp,1)==2
    max_u3 = max(abs(modeshape_disp(2,:)));
elseif size(modeshape_disp,1)==3
    max_u3 = max(abs(modeshape_disp(3,:)));
else
    error('size(modeshape,1) must equal to 2 or 3!');
end
max_x3 = max(abs(x3));
factor_scale = (max_x3/20)/max_u3*scale;

modeshape_u1 = modeshape_disp(1,:)*factor_scale;
% modeshape_u2 = modeshape_disp(2,:)*factor_scale;  % modeshape_u2 will not be used in this function
if size(modeshape_disp,1)==2
    modeshape_u3 = modeshape_disp(2,:)*factor_scale;
else
    modeshape_u3 = modeshape_disp(3,:)*factor_scale;
end

POS_x1 = zeros( length(x1), length(x3), length(Time) );
POS_x3 = zeros( length(x1), length(x3), length(Time) );
for t = 1:length(Time)
    pos_x1 = zeros( length(x1), length(x3) );
    pos_x3 = zeros( length(x1), length(x3) );
    for i = 1:num_x1
        for j = 1:num_x3
            u1 = modeshape_u1(j)*real(exp( 1i*xi*(x1(i)-cp*Time(t)) ));
            u3 = modeshape_u3(j)*real(exp( 1i*xi*(x1(i)-cp*Time(t)) ));
            pos_x1(i,j) = x1(i) + u1;
            pos_x3(i,j) = x3(j) + u3;
        end
    end
    POS_x1(:,:,t) = pos_x1;
    POS_x3(:,:,t) = pos_x3;
end
min_POS_x1 = min(min(min(POS_x1))); max_POS_x1 = max(max(max(POS_x1)));
min_POS_x3 = min(min(min(POS_x3))); max_POS_x3 = max(max(max(POS_x3)));

hf = figure('windowstyle','docked', 'color','w'); hold on;
ha = findobj(hf, 'type','axes');
if strcmp(nondeformedgrids, 'ON')    % plot the undeformed grids
    for i = 1:num_x1    % plot vertical line
        plot(ha, x1(i)*ones(size(x3)), x3, '-', 'color',[0.8,0.8,0.8], 'linewidth',1 );
    end
    for j = 1:num_x3    % plot horizontal line
        plot(ha, x1, x3(j)*ones(size(x1)), '-', 'color',[0.8,0.8,0.8], 'linewidth',1 );
    end
end
xlim([min_POS_x1,max_POS_x1]); ylim([min_POS_x3,max_POS_x3]); axis off;
text( (min_POS_x1+max_POS_x1)/2, min_POS_x3, 'Propagation direction (x1)', 'HorizontalAlignment','center', 'VerticalAlignment','cap');
text( min_POS_x1, (min_POS_x3+max_POS_x3)/2, 'Thickness (x3)', 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'rotation',90);

hx1 = []; hx3 = [];
for t = 1:length(Time)
    delete(hx1); delete(hx3);
    
    % plot the deformed grids
    hx1 = zeros(1,num_x1); hx3 = zeros(1,num_x3);
    for i = 1:num_x1    % plot vertical line
        hx1(i) = plot(ha, POS_x1(i,:,t), POS_x3(i,:,t), 'r-', 'linewidth',1 );
    end
    for j = 1:num_x3    % plot horizontal line
        hx3(j) = plot(ha, POS_x1(:,j,t), POS_x3(:,j,t), 'r-', 'linewidth',1 );
    end
    drawnow;
    
    if strcmp(savegif, 'ON')  % save gif file
        if ~exist('n_frame','var'), n_frame = 0; end  % initialize the number of frames for making gif
        n_frame = n_frame+1; frame = getframe(hf); im = frame2im(frame); [A, map] = rgb2ind(im, 256);
        if n_frame == 1
            imwrite(A, map, filename, 'gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(A, map, filename, 'gif','WriteMode','append','DelayTime',0.1);
        end
    end
end

end