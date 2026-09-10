function [varargout] = meshplot2D(pts,tri,nodeData,opt)
% Visualize multiple data on a 2D mesh.
%
% If the input "nodeData" contains multiple fields (such as a matrix where each
% column is a domain-continuous coeff, or a page-array where each page is 
% a simplex-continuous coeff), all data will be plotted on the same mesh and 
% displayed within the same figure. It is possible to change the displayed data 
% via the command fig.UserData.setIndex(...), or by using the mouse scroll wheel 
% while hovering over the figure.
%  
% SEE ALSO: PATCH, MESHPLOT3D
arguments
	pts (2,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	nodeData (:,:,:) double {mustBeFinite}
	opt.cmapLimits (1,2) {validateattributes(opt.cmapLimits,{'string','numeric'},{})} = ["",""];
	opt.fixedCLims (1,1) logical = false;
	opt.wrapScroll (1,1) logical = false;
	opt.titleIndex (1,1) logical = true;
	opt.plot3Dsurf (1,1) logical = true;
	opt.edgeColor = [];
	opt.faceColor = [];
	opt.Colormap = [];
end
[Ndim,Npts] = size(pts);
[~,Ntri] = size(tri);
tri = tri(1:Ndim+1,:);

% check/prepare nodeData
%===========================================================
nodeData = squeeze(nodeData);
if isempty(nodeData)
	plotMeshOnly = true;
	isdisjoint = false;
	Nframe = 1;
elseif ndims(nodeData)==3
	error('ERROR: the given node data has invalid size');
elseif any(size(nodeData)==Npts)
	plotMeshOnly = false;
	isdisjoint = false;
	if ((size(nodeData,1)~=Npts) & (size(nodeData,2)==Npts))
		nodeData = nodeData.';
	end
	Nframe = size(nodeData,2);
elseif any(size(nodeData)==Ntri)
	plotMeshOnly = false;
	isdisjoint = true;
	if (size(nodeData,1)~=Ntri) & (size(nodeData,2)==Ntri)
		nodeData = nodeData.';
	end
	Nframe = 1;
else
	error('ERROR: the given node data has invalid size');
end

% handle complex inputs
%-----------------------------------------------------------
isComplex = ~isreal(nodeData);
if isComplex
	colorData = angle(nodeData); % phase information
	nodeData = abs(nodeData).^2; % magnitude information
	opt.Colormap = hsv;
	opt.cmapLimits = (pi/2)*[-1,+1];
else
	colorData = nodeData;
end

% setup colormap/limits
%===========================================================
if isequal(opt.cmapLimits,["",""])
	if any(sign(nodeData(:))==+1) && any(sign(nodeData(:))==-1)
		opt.cmapLimits = ["avg-2*std","avg+2*std"];
	else
		opt.cmapLimits = ["min","max"];
	end
end

if isnumeric(opt.cmapLimits)
	cLims = opt.cmapLimits + zeros(Nframe,2);
elseif isstring(opt.cmapLimits)
	% choose between fixed and dynamic colormap limits
	cData = nodeData;
	if isdisjoint; cData = reshape(cData,[],size(cData,3)); end
	if opt.fixedCLims; cData = cData(:); end

	% get statistics for nodeData
	dmin = min(cData,[],1);
	dmax = max(cData,[],1);
	dmag = max(abs(dmin),abs(dmax));
	davg = mean(cData,1);
	dstd = std(cData)+eps;
	
	% evaluate input
	evalClim = @(x) feval(str2func("@(min,max,mag,avg,std)"+x),dmin,dmax,dmag,davg,dstd);
	cmin = evalClim(opt.cmapLimits(1));
	cmax = evalClim(opt.cmapLimits(2));
	cLims = [zeros(Nframe,1)+cmin(:),zeros(Nframe,1)+cmax(:)];
end
% make sure that clim(1) < clim(2)
issame = (cLims(:,1)==cLims(:,2));
isdecr = ~(cLims(:,1)<=cLims(:,2));
cLims(issame,:) = cLims(issame,:) + eps(cLims(issame,:)).*[-1,+1];
cLims(isdecr,:) = [cLims(isdecr,2),cLims(isdecr,1)];

% set colormap if unspecified
if isempty(opt.Colormap)
	if any(sign(cLims(:,1))==-1) & any(sign(cLims(:,2))==+1)
		% use divergent colormap
		opt.Colormap = cmap_dark;
	else
		% use monotonic colormap
		opt.Colormap = cmap_mono;
	end
end

% figure configs
%===========================================================
fig = figure('WindowScrollWheelFcn',@scroll_update);
fig.WindowStyle ='normal';
fig.Visible = 'off';
fig.Colormap = opt.Colormap;

% axis settings
ax = axes('Visible','off');
ax.NextPlot = 'replaceChildren';
ax.Units = 'normalized';
ax.InnerPosition = [0.04 0.01 0.92 0.92];

% set axes limits
xbnds = [min(pts(1,:)), max(pts(1,:))];
ybnds = [min(pts(2,:)), max(pts(2,:))];
xybox = max([diff(xbnds),diff(ybnds)]) * [-1/2, 1/2];
ax.XLim = mean(xbnds)+xybox;
ax.YLim = mean(ybnds)+xybox;
ax.Colormap = opt.Colormap;

% set labeling method
%-----------------------------------------------------------
if opt.titleIndex
	% use whole indices
	label = @(x) ['   ','[',sprintf('%i',x),']'];
else
	% use fractional indices
	label = @(x) ['   ','[',sprintf('%.3f',(x-1)/(Nframe-1)),']'];
end

% setup displayed data index
%-----------------------------------------------------------
indx = 1;
indxData.setIndex = @update_axis;
set(fig,'UserData',indxData)

% populate frames
%===========================================================
if plotMeshOnly
	% plot boundary mesh and exit
	if isempty(opt.faceColor); opt.faceColor = '#9a8566'; end
	if isempty(opt.edgeColor); opt.edgeColor = '#000000'; end

	patch('Vertices',pts.','Faces',tri.','Parent',ax, ...
		'FaceColor',opt.faceColor,'FaceAlpha',1,'EdgeColor',opt.edgeColor,'EdgeAlpha',1)
	set(fig,"Visible","on");
	if nargout > 0
		varargout{1} = fig;
	end
	return
else
	if isempty(opt.edgeColor); opt.edgeColor = '#242424'; end
end

if isdisjoint
	% setup axes for disjointed data
	zpts = [];
	ztri = zeros(Ndim+1,Ntri);
	for ss = 1:Ntri
		ztri(:,ss) = size(zpts,2) + [1:Ndim+1].';
		zpts = [zpts,[pts(:,tri(:,ss)); nodeData(ss,:)]];
	end
	ax.PlotBoxAspectRatio = [1 1 0.5];
	ax.ZLim = [min(cLims(:,1)),max(cLims(:,2))];

	pat = patch('Vertices',zpts.','Faces',ztri.','Parent',ax, ...
		'EdgeColor',opt.edgeColor,'FaceColor','interp',"FaceVertexCData",zpts(3,:).');
elseif opt.plot3Dsurf
	% setup axes for 3D surface display
	vertX = @(ind) [pts.',nodeData(:,ind)];
	ax.PlotBoxAspectRatio = [1 1 0.5];
	ax.ZLim = [min(cLims(:,1)),max(cLims(:,2))];

	pat = patch('Vertices',vertX(indx),'Faces',tri.','Parent',ax, ...
		'FaceColor','interp','EdgeColor',opt.edgeColor);
	set(pat,'Cdata',colorData(:,indx))
else
	% setup axes for 2D surface display
	ax.PlotBoxAspectRatio = [1 1 1];
	pat = patch('Vertices',pts.','Faces',tri.','Parent',ax, ...
		'FaceColor','interp','EdgeColor',opt.edgeColor);
	if isComplex
		set(pat,'AlphaDataMapping','scaled');
		set(pat,'FaceVertexAlphaData',nodeData(:,indx));
	end
	set(pat,'Cdata',colorData(:,indx))
end

% display figure and define update functions
%-----------------------------------------------------------
set(fig,'Name',label(indx)); % set figure title
set(ax,'CLim',cLims(indx,:)); % set figure clims
set(fig,"Visible","on");
if nargout > 0
	varargout{1} = fig;
end

%-----------------------------------------------------------
	function scroll_update(~,event)
	% use mouse scroll wheel to update the displayed data index
	new_indx = indx + (event.VerticalScrollCount);
	update_axis(new_indx)
	end
%-----------------------------------------------------------
	function update_axis(new_indx)
	% updates the displayed data after performing checks
	doUpdate = checkIndex(new_indx);

	if ~doUpdate; return; end
	set(pat,'Cdata',colorData(:,indx));

	if opt.plot3Dsurf
		set(pat,'Vertices',vertX(indx));
	elseif isComplex
		set(pat,'FaceVertexAlphaData',nodeData(:,indx));
	end

	set(fig,'Name',label(indx)); % update figure title
	set(ax,'CLim',cLims(indx,:)); % update color-data limits
	end
%-----------------------------------------------------------
	function [doUpdate] = checkIndex(new_indx)
	% validates new_indx, enforces indexing bounds, and determines if update is needed
	doUpdate = false;
	if isnumeric(new_indx) && (round(new_indx)==new_indx)
		if opt.wrapScroll
			new_indx = mod(new_indx, Nframe);
			if new_indx==0; new_indx = Nframe; end
		else
			new_indx = min(new_indx,Nframe);
			new_indx = max(new_indx,1);
		end
		if new_indx~=indx
			indx = new_indx;
			doUpdate = true;
		end
	end
	end
%-----------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%