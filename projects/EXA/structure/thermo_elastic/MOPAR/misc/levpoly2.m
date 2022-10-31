function [levelset,ex,ey,nx,ny,d_distf]=levpoly2(func,nx,ny,dc,pflag,s,ds)
%function [levelset,ex,ey,nx,ny,d_distf]=levpoly(func,nx,ny,dc,pflag,s,ds)
%--------------------------------------------------------------------------
%
%  generates signed distance levelset field from ploygon
%
%  input: func  - function for generating splines
%                 call [cord,topo,prop,numnod,numele]=func(s)
%                 cord -  numnod x 2 (x,y)
%                 topo -  numele x 3 (edge type,node1, node2)
%                 prop -  numele x 1 (polygon id)
%         nx,ny - mesh dimensions (cells - not nodes)
%         dc    - size of square cell
%         pflag - plot flag
%         s     - vector of optimization varialbes
%                 need to fit to design file
%                 size: numvar x 1
%         ds    - matrix of derivatives
%                 size: numvar x numder
%
%  output: 
%--------------------------------------------------------------------------

colvec=['b','g','c','m','y','r'];

seps=1e-5;

% generate mesh with lower left corner at (0,0)

tic

[ex,ey] = meshgrid(0:nx,0:ny);
ex      = dc*ex;
ey      = dc*ey;
ez      = 0*ex;

fprintf('\n\n mpoly - generated mesh : %e sec\n',toc);

% generate bondary polygon of obstacle

tic

[cord,topo,prop,numnod,numele]=feval(func,s);

numpoly=max(prop);

% derivative of boundary nodes by finite difference

numder=size(ds,2);
if numder >= 1
    numder=1;
end

d_cord=zeros(numder,numnod,2);

for ider=1:numder
    sp=s+seps*ds(:,ider);
    sm=s-seps*ds(:,ider);
    [cordp]=feval(func,sp);
    [cordm]=feval(func,sm);
    d_cord(ider,:,:)=(cordp-cordm)/2/seps;
end

% plot boundary

if pflag
    figure(1)
    for ie=1:numele
        inod=topo(ie,2:3);
        cl=[colvec(prop(ie)) '-'];
        plot(cord(inod,1),cord(inod,2),cl,'LineWidth',1);
        hold on
    end
    mesh(ex,ey,ez);
    axis equal
end

fprintf('\n\n mpoly - generated splines : %e sec\n',toc);

% evaluate distance field

tic

distf=1e9*ones(nx+1,ny+1);
% dispol=zeros(nx+1,ny+1);
d_distf=zeros(numder,nx+1,ny+1);

kp=0;

for ip=1:numpoly
    
    numelp=length(find(prop == ip));
    pvec=zeros(2*numelp,2);

    for ie=1:numelp
        pvec(2*ie-1,:)=cord(topo(kp+ie,2),:);
        pvec(2*ie  ,:)=cord(topo(kp+ie,3),:);
    end

    iact=zeros(nx+1,ny+1);
    xp=reshape(ex',[],1);
    yp=reshape(ey',[],1);
    
    [xy,distance] = distance2curve(cord,[xp yp]);
    inp=inpolygon(xp, yp, cord(:,1),cord(:,2));
    distf=reshape((1-2*inp).*distance,nx+1,ny+1);
    iact(:,:)=1;
    
%     for ix=1:nx+1
%         for iy=1:ny+1
%             dsf=p_poly_dist(dc*ix-dc,dc*iy-dc,pvec(:,1),pvec(:,2));
%             if abs(dsf) < abs(distf(ix,iy))
%                 distf(ix,iy)=dsf;
%                 dispol(ix,iy)=ip;
%                 iact(ix,iy)=1;
%             end
%         end
%     end
    
    pvecp=zeros(2*numelp,2);
    pvecm=zeros(2*numelp,2);

    for ider=1:numder    
        for ie=1:numelp
            pvecp(2*ie-1,:)=cord(topo(kp+ie,2),:)+seps*squeeze(d_cord(ider,topo(kp+ie,2),:))';
            pvecp(2*ie  ,:)=cord(topo(kp+ie,3),:)+seps*squeeze(d_cord(ider,topo(kp+ie,3),:))';
            pvecm(2*ie-1,:)=cord(topo(kp+ie,2),:)-seps*squeeze(d_cord(ider,topo(kp+ie,2),:))';
            pvecm(2*ie  ,:)=cord(topo(kp+ie,3),:)-seps*squeeze(d_cord(ider,topo(kp+ie,3),:))';
        end

        for ix=1:nx+1
            for iy=1:ny+1
                if iact(ix,iy) > 0
                    distfp=p_poly_dist(dc*ix-dc,dc*iy-dc,pvecp(:,1),pvecp(:,2));
                    distfm=p_poly_dist(dc*ix-dc,dc*iy-dc,pvecm(:,1),pvecm(:,2));
                    d_distf(ider,ix,iy) = (distfp-distfm)/2/seps;
                end
            end
        end
    end
    
    kp=kp+numelp;
end
    
if pflag
    figure(2);
    xm     = zeros(nx+1,ny+1);
    ym     = zeros(nx+1,ny+1);
    distfp = zeros(nx+1,ny+1);
    for ix=1:nx+1
        for iy=1:ny+1
          xm(ix,iy)=dc*ix-dc;
          ym(ix,iy)=dc*iy-dc;
          distfp(ix,iy)=distf(ix,iy);
        end
    end
    for ie=1:numele
        inod=topo(ie,2:3);
        plot(cord(inod,1),cord(inod,2),'k-');
        hold on
    end
    contour(xm,ym,distfp);
    hold on
    contour(xm,ym,distfp,[0 0],'r','LineWidth',2);
    colorbar
    axis equal
    
%     for ix=1:nx+1
%         for iy=1:ny+1
%           xp=dc*ix-dc;
%           yp=dc*iy-dc;
%           plot(xp,yp,'k.');
%           txt=sprintf('(%d,%d):%f',ix,iy,distf(ix,iy));
%           text(xp,yp,txt,'FontSize',8);
%         end
%     end

end

fprintf('\n mpoly - built distance field : %e sec\n',toc);

levelset=reshape(distf,(nx+1)*(ny+1),1);
%%% save levelset.mat levelset nx ny ex ey;

