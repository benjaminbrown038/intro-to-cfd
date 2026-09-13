function result = cfd_stokes(edges, mode, speed, mu, obstacle, penalty)
% Staggered finite-volume steady Stokes on bounded 2D/3D tensor grids.
% Same equations and boundary rules as common/python/stokes.py.
if nargin<5, obstacle=false; end
if nargin<6, penalty=1e4; end
start=tic; dim=numel(edges);
assert(any(dim==[2,3]) && mu>0 && isfinite(mu));
assert(strcmp(mode,'cavity') || (strcmp(mode,'channel') && dim==2));
n=zeros(1,dim); h=cell(1,dim); centers=h;
for d=1:dim
    edges{d}=edges{d}(:)';
    assert(numel(edges{d})>=4 && all(isfinite(edges{d})) && all(diff(edges{d})>0));
    n(d)=numel(edges{d})-1; h{d}=diff(edges{d});
    centers{d}=(edges{d}(2:end)+edges{d}(1:end-1))/2;
end
volumes=ones(n); widthsgrid=cell(1,dim); [widthsgrid{:}]=ndgrid(h{:});
for d=1:dim, volumes=volumes.*widthsgrid{d}; end
np=prod(n); continuity=zeros(np,1);
for row=1:np
    idx=unpack(n,row); point=zeros(1,dim);
    for d=1:dim, point(d)=centers{d}(idx(d)); end
    for d=1:dim
        area=volumes(row)/h{d}(idx(d));
        for side=0:1
            if (side==0 && idx(d)==1) || (side==1 && idx(d)==n(d))
                bp=point; bp(d)=edges{d}(1+side*n(d));
                continuity(row)=continuity(row)-(2*side-1)*area*wall(d,d,side,bp,mode,speed,edges);
            end
        end
    end
end
assert(abs(sum(continuity))<1e-10,'Boundary flux mismatch');
blocks=cell(1,dim); loads=blocks; shapes=blocks; masks=blocks; dualvolumes=blocks;
br=[]; bc=[]; bv=[]; offset=0;
for comp=1:dim
    shape=n; shape(comp)=shape(comp)-1; shapes{comp}=shape;
    coords=centers; coords{comp}=edges{comp}(2:end-1);
    widths=h; widths{comp}=(h{comp}(2:end)+h{comp}(1:end-1))/2;
    nr=prod(shape); rr=zeros(nr*(2*dim+1),1); cc=rr; vv=rr; nz=0;
    rhs=zeros(nr,1); mask=rhs; dualvol=rhs;
    for row=1:nr
        idx=unpack(shape,row); point=zeros(1,dim); vol=1;
        for d=1:dim, point(d)=coords{d}(idx(d)); vol=vol*widths{d}(idx(d)); end
        dualvol(row)=vol; diagonal=0;
        for d=1:dim
            area=vol/widths{d}(idx(d));
            for side=0:1
                neighbor=idx; neighbor(d)=neighbor(d)+2*side-1;
                if neighbor(d)>=1 && neighbor(d)<=shape(d)
                    distance=abs(coords{d}(neighbor(d))-point(d));
                    conductance=mu*area/distance;
                    nz=nz+1; rr(nz)=row; cc(nz)=pack(shape,neighbor); vv(nz)=-conductance;
                else
                    bound=edges{d}(1+side*n(d)); distance=abs(point(d)-bound);
                    conductance=mu*area/distance; bp=point; bp(d)=bound;
                    rhs(row)=rhs(row)+conductance*wall(comp,d,side,bp,mode,speed,edges);
                end
                diagonal=diagonal+conductance;
            end
        end
        if strcmp(mode,'channel') && comp==1
            height=edges{2}(end)-edges{2}(1); rhs(row)=rhs(row)+8*mu*speed/height^2*vol;
        end
        if obstacle
            radius=.15*(edges{2}(end)-edges{2}(1));
            cx=(edges{1}(end)+edges{1}(1))/2; cy=(edges{2}(end)+edges{2}(1))/2;
            mask(row)=double((point(1)-cx)^2+(point(2)-cy)^2<=radius^2);
            diagonal=diagonal+penalty*vol*mask(row);
        end
        nz=nz+1; rr(nz)=row; cc(nz)=row; vv(nz)=diagonal;
        left=idx; right=idx; right(comp)=right(comp)+1; area=1;
        for d=1:dim, if d~=comp, area=area*h{d}(idx(d)); end, end
        br(end+1:end+2,1)=[pack(n,left);pack(n,right)];
        bc(end+1:end+2,1)=[offset+row;offset+row]; bv(end+1:end+2,1)=[area;-area];
    end
    blocks{comp}=sparse(rr(1:nz),cc(1:nz),vv(1:nz),nr,nr);
    loads{comp}=rhs; masks{comp}=mask; dualvolumes{comp}=dualvol; offset=offset+nr;
end
K=blkdiag(blocks{:}); B=sparse(br,bc,bv,np,offset); Bt=B(2:end,:);
A=[K,-Bt';-Bt,sparse(np-1,np-1)]; rhs=[vertcat(loads{:});-continuity(2:end)];
assembly_seconds=toc(start); timer=tic; sol=A\rhs; solve_seconds=toc(timer);
assert(all(isfinite(sol)),'Nonfinite solution');
velocity=sol(1:offset); pressure=reshape([0;sol(offset+1:end)],n);
divergence=reshape((B*velocity-continuity)./volumes(:),n);
faces=cell(1,dim); cell_velocity=faces; pos=0; force=zeros(1,dim);
for comp=1:dim
    shape=shapes{comp}; count=prod(shape); values=velocity(pos+1:pos+count); pos=pos+count;
    fs=n; fs(comp)=fs(comp)+1; full=zeros(fs);
    sl=repmat({':'},1,dim); sl{comp}=2:n(comp); full(sl{:})=reshape(values,shape);
    % Set prescribed normal velocities on boundary faces.
    for row=1:numel(full)
        idx=unpack(fs,row);
        if idx(comp)==1 || idx(comp)==fs(comp)
            side=double(idx(comp)==fs(comp)); point=zeros(1,dim);
            for d=1:dim
                if d==comp, point(d)=edges{d}(1+side*n(d)); else, point(d)=centers{d}(idx(d)); end
            end
            full(row)=wall(comp,comp,side,point,mode,speed,edges);
        end
    end
    a=repmat({':'},1,dim); b=a; a{comp}=1:n(comp); b{comp}=2:n(comp)+1;
    cell_velocity{comp}=.5*(full(a{:})+full(b{:})); faces{comp}=full;
    if obstacle, force(comp)=sum(penalty*masks{comp}.*dualvolumes{comp}.*values); end
end
metrics=struct('cells',np,'unknowns',size(A,1),'assembly_seconds',assembly_seconds,'solve_seconds',solve_seconds,...
    'relative_residual',norm(A*sol-rhs)/max(norm(rhs),1e-30),'max_divergence',max(abs(divergence(:))),'penalty_force',force);
result=struct('edges',{edges},'centers',{centers},'velocity',{cell_velocity},'faces',{faces},'pressure',pressure,...
    'divergence',divergence,'volumes',volumes,'metrics',metrics);
end

function idx=unpack(shape,row)
c=cell(1,numel(shape)); [c{:}]=ind2sub(shape,row); idx=cell2mat(c);
end
function row=pack(shape,idx)
c=num2cell(idx); row=sub2ind(shape,c{:});
end
function value=wall(comp,axis,side,point,mode,speed,edges)
value=0;
if strcmp(mode,'cavity') && comp==1 && axis==2 && side==1
    value=speed;
    for d=1:numel(edges)
        if d~=2, value=value*sin(pi*(point(d)-edges{d}(1))/(edges{d}(end)-edges{d}(1)))^2; end
    end
elseif strcmp(mode,'channel') && comp==1 && axis==1
    y=(point(2)-edges{2}(1))/(edges{2}(end)-edges{2}(1)); value=4*speed*y*(1-y);
end
end
