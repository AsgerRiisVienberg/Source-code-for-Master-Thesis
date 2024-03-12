function enhance_plot(fontname,fontsize,linewid,markersiz,includeLegend)

%  enhance_plot([fontname,fontsize]);
%
%  Function to enhance MATLAB's lousy text choices on plots.  Sets the
%  current figure's Xlabel, Ylabel, Title, and all Text on plots, plus
%  the axes-labels to the "fontname" and "fontsize" input here where
%  the defaults have been set to 'times' and 16.
%  Also sets all plotted lines to "linewid" and all markers to size
%  "markersiz".  The defaults are 2 and 8.
%
%  INPUTS:  fontname:   (Optional,DEF='TIMES') FontName string to use
%                       MATLAB's ugly default is 'Helvetica'
%           fontsize:   (Optional,DEF=16) FontSize integer to use
%                       MATLAB's tiny default is 10
%           linewid:    (Optional,DEF=2) LineWidth integer to use
%                       MATLAB's skinny default is 0.5
%           markersiz:  (Optional,DEF=8) MarkerSize integer to use
%                       MATLAB's squinty default is 6
%  for all inputs, if pass 0, use default
%                 if pass -1, use MATLAB's default
%
% Modifications
%  6-Sep-2022 A. Vienberg
%       updated the function to be compatible with MATLAB 2022a
% 
%  1-Nov-2015 M. Jerez
%       added latex interpreter and automatic correction of legends
%
%  19-Feb-2002 J. Nelson
%       added linewid and markersiz to help squinting readers
%
%  20-Feb-2002 J. Nelson
%       added check for legend.  If legend exists, increase the 
%           line and marker size, also increase the font to 
%           fontsize-2 (2 points smaller than title and labels)
%  25-Feb-2002 J. Nelson
%       added lgd (legend) input check to fix legend problems.
%======================================================================




if (~exist('fontname','var')||(fontname==0))
  fontname = 'times';
elseif (fontname==-1)
  fontname = 'helvetica';
end
if (~exist('fontsize','var')||(fontsize==0))
  fontsize = 16;
elseif (fontsize==-1)
  fontsize=10;
end
if (~exist('linewid','var')||(linewid==0))
  linewid=2;
elseif (linewid==-1)
  linewid=0.5;
end
if (~exist('markersiz','var')||(markersiz==0))
  markersiz = 8;
elseif (markersiz==-1)
  markersiz = 6;
end
if (~exist('legendSetting','var')||(includeLegend==0))
    includeLegend = true;
elseif (includeLegend == -1)
    includeLegend = false;
end

if (~exist('lgd')|(lgd==0)|(lgd==-1))
 lgd=0;
end

Hf=gcf;
Ha=gca;
Hx=get(Ha,'XLabel');
Hy=get(Ha,'YLabel');
Ht=get(Ha,'Title');
set(Ha,'LineWidth',.75);
set(Hx,'fontname',fontname);
set(Hx,'fontsize',fontsize);
set(Hx,'Interpreter','latex');
set(Hy,'fontname',fontname);
set(Hy,'fontsize',fontsize);
set(Hy,'Interpreter','latex');
set(Ha,'fontname',fontname);
set(Ha,'fontsize',fontsize);

set(Ha,'YaxisLocation','right')
set(Ha,'YaxisLocation','left')
set(Ht,'fontname',fontname);
set(Ht,'fontsize',fontsize);
set(Ht,'Interpreter','latex');
set(Hy,'VerticalAlignment','bottom');
set(Hx,'VerticalAlignment','cap');
set(Ht,'VerticalAlignment','baseline');
Hn = get(Ha,'Children');
n = length(Hn);
if n > 0
  typ = get(Hn,'Type');
  for j = 1:n
    if strcmp('text',typ(j,:))
      set(Hn(j),'fontname',fontname);
      set(Hn(j),'fontsize',fontsize);
      set(Hn(j),'Interpreter','latex');

    end
    if strcmp('line',typ(j,:))
      set(Hn(j),'LineWidth',linewid);
      set(Hn(j),'MarkerSize',markersiz);
    end
  end
end
% 
if lgd ~= 0
    legh=legend;
    Hn=get(legh,'String');
    n = length(Hn);
    if n > 0
        set(legh,'fontname',fontname);
        set(legh,'fontsize',fontsize-2);
        set(legh,'Interpreter','latex');
    end
end
figure(Hf);