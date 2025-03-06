function [argval,outputs]=parseargpair(inputs,name,default)

if(nargin<3 || ~exist('default','var'))
  default=[];
end

argnum=find(strcmp(name,inputs));
if(~isempty(argnum))
  argval=inputs{argnum(1)+1};
  inputs(argnum(1)+(0:1))=[];
else
  argval=default;
end

outputs=inputs;

end