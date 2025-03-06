% function ff=newfig(name,[param-pairs])
% 
% creates a new figure with default publication parameter
% 
% if you give a specific name to the figure, next time you call the
% function the existing figures with the same name will be first closed.
%
% ['param',value] lets you set param-value pairs to the figure according to
%                   the figure properties

function [ff,aa]=newfig(name,varargin)

publicationparam={
    'DefaultAxesTickDirMode'        'manual'
    'DefaultAxesTickDir'            'out'
    'DefaultLineLineWidth'          1
    'DefaultAxesLineWidth'          1
    'DefaultAxesTickLength'         [0.015 0.015]
    'DefaultAxesFontSize'           8
    'DefaultTextFontSize'           8
    'DefaultAxesTitleFontSize'      1.2
    'DefaultAxesLabelFontSize'      1.2
    'DefaultLineMarkerSize'         15
    'Unit'             'centimeter'
    'Position'         [17 14 8.5 6.3]
    }';

if(~exist('name','var') || isempty(name))
    name=num2str(round(rand*1e6));
else
    cprintf('_comments','%s\n',name)
    old=findobj('newfig_name',name);
    delete(old);
end

[subplt,varargin]=parseargpair(varargin,'subplot',[]);

f=figure('name',name,publicationparam{:},varargin{:});
if(~isempty(subplt))
  ax=subplot_ax(subplt);
end

hProp = addprop(f,'newfig_name');
set(f,'newfig_name',name);
hProp.SetAccess = 'private';

if(nargout>0)
    ff=f;
    if(~isempty(subplt))
      aa=ax;
    end
end