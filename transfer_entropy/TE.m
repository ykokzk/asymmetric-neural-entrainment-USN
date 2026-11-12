%{ 
This script is a class function for calculating transfer entropy between
input variables.

Author: H.Yokoyama

input: 
    x, y, xdelay, ydelay: row vector with size [1 x N]  

example of use:
    >> x = rand(1,1001);
    >> y = randn(1,1001);
    >>
    >> X = x(1:end-1);
    >> Y = y(1:end-1);
    >> Xdelay = x(2:end);
    >> Ydelay = y(2:end);
    >>
    >> te = TE(X,Y,Xdelay, Ydelay);
    >> te.get_PDF()
    >> [TEy2x, TEx2y] = te.get_TE()
%}

classdef TE < handle
    properties
        x;          % input variable x
        y;          % input variable y
        xdelay;     % input variable x with delay
        ydelay;     % input variable y with delay
        Nbins;      % optimal bin number 
        px;         % probability density function p(x)
        py;         % probability density function p(y)
        pxy;        % joint probability density function p(x, y)
        pyx;        % joint probability density function p(y, x)
        pxx;        % joint probability density function p(xdelay, x)
        pyy;        % joint probability density function p(ydelay, y)
        pxxy;       % joint probability density function p(xdelay, x, y)
        pyyx;       % joint probability density function p(ydelay, y, x)
        TEy2x;      % transfer entropy from y to x
        TEx2y;      % transfer entropy from x to y
    end

    
    methods
        function obj = TE(x, y, xdelay, ydelay)
            obj.x      = x;
            obj.y      = y;
            obj.xdelay = xdelay;
            obj.ydelay = ydelay;
            
            N          = length(x);
            Nbins      = 1 + 3.322*log(N);
            obj.Nbins  = round(Nbins);

        end

        %%
        function get_PDF(obj)
        % Wiritten by Yokoyama
        %
            Bininfo.x      = [obj.Nbins; min(obj.x); max(obj.x)];
            Bininfo.y      = [obj.Nbins; min(obj.y); max(obj.y)];
            Bininfo.xdelay = [obj.Nbins; min(obj.xdelay); max(obj.xdelay)];
            Bininfo.ydelay = [obj.Nbins; min(obj.ydelay); max(obj.ydelay)];
        
            obj.px     = histndim(obj.x.',                                Bininfo.x,                                        "p"); %Nbinsで定義された各ビンにおいて、各ビンに入るデータｘの確率を返す?xの確率分布
            obj.py     = histndim(obj.y.',                                Bininfo.y,                                        "p"); 
            obj.pxy   = histndim([obj.x.', obj.y.'],                   [Bininfo.x, Bininfo.y],                        "p"); %x,yの結合確率分布
            obj.pyx   = histndim([obj.y.', obj.x.'],                   [Bininfo.y, Bininfo.x],                        "p");
            obj.pxx   = histndim([obj.xdelay.', obj.x.'],           [Bininfo.xdelay, Bininfo.x],                "p"); 
            obj.pyy   = histndim([obj.ydelay.', obj.y.'],           [Bininfo.ydelay, Bininfo.y],                "p");
            obj.pxxy = histndim([obj.xdelay.', obj.x.', obj.y.'], [Bininfo.xdelay, Bininfo.x, Bininfo.y], "p");
            obj.pyyx = histndim([obj.ydelay.', obj.y.', obj.x.'], [Bininfo.ydelay, Bininfo.y, Bininfo.x], "p");
        end
        
        function [TEy2x, TEx2y, Hx, Hxx, Hxy, Hxxy] = get_TE(obj)
            Hx    = double(-obj.px(obj.px~=0).'    * log2(obj.px(obj.px~=0)));
            Hxx   = double(-obj.pxx(obj.pxx~=0).'  * log2(obj.pxx(obj.pxx~=0)));
            Hxy   = double(-obj.pxy(obj.pxy~=0).'  * log2(obj.pxy(obj.pxy~=0)));
            Hxxy  = double(-obj.pxxy(obj.pxxy~=0).'* log2(obj.pxxy(obj.pxxy~=0)));
            
            TEy2x = (Hxx + Hxy - Hx - Hxxy);

            if TEy2x < 0 || isnan(TEy2x) || abs(TEy2x)==Inf
                TEy2x = 0;
            end
            
            Hy    = double(-obj.py(obj.py~=0).'    * log2(obj.py(obj.py~=0)));
            Hyy   = double(-obj.pyy(obj.pyy~=0).'  * log2(obj.pyy(obj.pyy~=0)));
            Hyx   = double(-obj.pyx(obj.pyx~=0).'  * log2(obj.pyx(obj.pyx~=0)));
            Hyyx  = double(-obj.pyyx(obj.pyyx~=0).'* log2(obj.pyyx(obj.pyyx~=0)));

            TEx2y = (Hyy + Hyx - Hy - Hyyx);
            
            if TEx2y < 0 || isnan(TEx2y) || abs(TEx2y)==Inf
                TEx2y = 0;
            end
            obj.TEy2x = TEy2x;
            obj.TEx2y = TEx2y;
        end
    end
    %%%%%%%%%%%%%%

end

function [v,t] = histndim(x,b,mode)
%HISTNDIM - generates and/or plots an n-dimensional histogram
%
%  Inputs:  X(m,d)   is the input data: each row is one d-dimensiona data point
%           B(3,d)   specifies the histogram bins.
%                         B(1,:) gives the number of bins in each dmension [default 10]
%                         B(2,:) gives the minimum of the first bin in each dimension [default min(X)]
%                         B(3,:) gives the maximum of the last bin in each dimension [default max(X)]
%                    If B has only one column, the same values are use for al dimensions 
%                    If B(1,i)=0 then that dimension will be ignored (and excluded from V)
%           MODE     is a character string containing a combination of the following:
%                        'z' for zero base in the 2D plot [default base = min(V)]
%                        'p' to scale V as probabilities [default actual counts]
%                        'h' to plot a histogam even if output arguments are present
%
% Outputs:  V        d-dimensional array containing the histogram counts
%           T        d-element cell array. d{i} contains the bin boundary values for
%                    the i'th dimension. The length of d{i} is one more than the number of bins
%                    in that dimension.
%
%                    Note that if any of B(1,:) are zero then the number of dimensions in V and elements
%                    of T will be correspondingly reduced.
%
% Example: histndim(randn(100000,2),[20 -3 3]','pz');

%	   Copyright (C) Mike Brookes 2004
%      Version: $Id: histndim.m 713 2011-10-16 14:45:43Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [n,d]=size(x);
    if nargin<3
        mode=' ';
        if(nargin<2)
            b=repmat(10,1,d);
        end
    end

    if size(b,2)==1
        b=repmat(b,1,d);
    end
    if size(b,1)<3
        mi=min(x,[],1);
        ma=max(x,[],1);
        w=(ma-mi)./(b(1,:)-0.001);  % nudge slightly to make sure al points included
        b(3,:)=ma+0.0005*w;
        b(2,:)=mi-0.0005*w;
    end

    acd=find(b(1,:)>0);
    sv=b(1,acd);
    nbt=prod(sv);
    t=cell(length(acd),1);

    % loop through each dimension
    k=1;        % indexing factor
    ok=repmat(1>0,n,1);
    ix=repmat(nbt-sum(cumprod(sv)),n,1);
    for i=1:length(acd)
        j=acd(i);
        bw=b(1,j)/(b(3,j)-b(2,j));
        bi=ceil((x(:,j)-b(2,j))*bw);
        ok=ok & (bi>0) & (bi<=b(1,j));
        ix(ok)=ix(ok)+k*bi(ok);
        k=k*b(1,j);
        t{i}=b(2,j)+(0:b(1,j))/bw;
    end
    v=full(sparse(ix(ok),1,1,nbt,1));
    if length(sv)>1
        v=reshape(v,sv);
    end
    if any(mode=='p')
        v=v/n;
    end

    if ~nargout || any(mode=='h')
        svg=find(sv>1);
        if length(svg)==1
            j=acd(svg);
            bar(b(2,j)+(0.5:sv(svg)-0.5)*(b(3,j)-b(2,j))/b(1,j),v(:));
        elseif length(svg)==2
            j=acd(svg(1));
            k=acd(svg(2));
            bj=b(1,j);
            bk=b(1,k);
            %     imagesc(b(2:3,k),b(2:3,j),reshape(v,b(1,j),b(1,k)));
            vda=kron(reshape(v,bj,bk),[1 1;1 1]);
            if any(mode=='z')
                ba=0;
            else
                ba=min(vda(:));
            end
            vda=[repmat(ba,1,2*bk+2);repmat(ba,2*bj,1) vda repmat(ba,2*bj,1);repmat(ba,1,2*bk+2)];
            jda=kron(t{svg(1)},[1 1]);
            jda=jda-(jda(3)-jda(2))*0.01*[-0.5 (-1).^(1:2*bj) 0.5]; % nudge slightly to avoid MATLAB plotting bug
            kda=kron(t{svg(2)},[1 1]);
            kda=kda-(kda(3)-kda(2))*0.01*[-0.5 (-1).^(1:2*bk) 0.5];
            surf(jda,kda,vda');
            ylabel(sprintf('Axis %d',k));
            xlabel(sprintf('Axis %d',j));
            colorbar;
        else
            fprintf(2,'Error in %s: Cannot plot 3-D histogram\n',mfilename);
        end
    end
end
    