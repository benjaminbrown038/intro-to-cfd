function cfd_run_cases(tier,base)
% Drivers for Tiers 3-5. Each tier's run_case.m adds this common directory.
names={'cavity','channel','adaptive'};
c=jsondecode(fileread(fullfile(base,'cases',[names{tier-2},'.json'])));
out=fullfile(base,'results','matlab'); if ~exist(out,'dir'), mkdir(out); end
if tier==3
    e=linspace(0,1,c.cells+1);
    r=cfd_stokes({e,e},'cavity',c.lid_speed,c.viscosity);
    save_result(r,out,'cavity_2d');disp(r.metrics);
elseif tier==4
    rows=struct([]);
    for k=1:numel(c.cell_counts)
        n=c.cell_counts(k); edges={linspace(0,2,2*n+1),linspace(0,1,n+1)};
        r=cfd_stokes(edges,'channel',c.peak_speed,c.viscosity);
        y=r.centers{2}; exact=4*c.peak_speed*y.*(1-y);
        err=r.velocity{1}-repmat(exact,numel(r.centers{1}),1);
        rows(k)=struct('cells',r.metrics.cells,'ny',n,'velocity_rms_error',sqrt(mean(err(:).^2)),'max_divergence',r.metrics.max_divergence);
        save_result(r,out,sprintf('channel_%d',n));
    end
    writetable(struct2table(rows),fullfile(out,'comparison.csv'));
    obstacle=cfd_stokes(edges,'channel',c.peak_speed,c.viscosity,true,c.penalty);
    save_result(obstacle,out,'penalized_cylinder');disp(struct2table(rows));
else
    rows=struct([]); q=0;
    for dim=2:3
        ref=cfd_stokes(repmat({linspace(0,1,c.reference_cells+1)},1,dim),'cavity',c.lid_speed,c.viscosity);
        save_result(ref,out,sprintf('reference_%dd',dim));
        edges=repmat({linspace(0,1,c.initial_cells+1)},1,dim);
        for cycle=0:c.cycles-1
            a=cfd_stokes(edges,'cavity',c.lid_speed,c.viscosity);
            uniform=cellfun(@(e)linspace(e(1),e(end),numel(e)),edges,'UniformOutput',false);
            u=cfd_stokes(uniform,'cavity',c.lid_speed,c.viscosity);
            save_result(a,out,sprintf('adaptive_%dd_cycle_%d',dim,cycle));
            save_result(u,out,sprintf('uniform_%dd_cycle_%d',dim,cycle));
            both={a,u}; types={'adaptive','uniform'};
            for j=1:2
                r=both{j};q=q+1;
                rows(q)=struct('dimension',dim,'cycle',cycle,'mesh',types{j},'cells',r.metrics.cells,'unknowns',r.metrics.unknowns,...
                    'sample_relative_error',sample_error(r,ref),'solve_seconds',r.metrics.solve_seconds,'max_divergence',r.metrics.max_divergence);
            end
            edges=refine_grid(a,c.add_per_axis);
        end
    end
    writetable(struct2table(rows),fullfile(out,'comparison.csv'));
    f=figure('Visible','off');
    for dim=2:3
        subplot(1,2,dim-1);hold on;
        for name={'adaptive','uniform'}
            selected=rows([rows.dimension]==dim & strcmp({rows.mesh},name{1}));
            plot([selected.unknowns],[selected.sample_relative_error],'o-','DisplayName',name{1});
        end
        xlabel('Unknowns');ylabel('Relative interior sample error');title(sprintf('%dD',dim));legend('show');grid on;
    end
    print(f,fullfile(out,'adaptive_comparison.png'),'-dpng','-r150');close(f);disp(struct2table(rows));
end
save(fullfile(out,'case.mat'),'c');
end

function save_result(r,out,name)
save(fullfile(out,[name,'.mat']),'r');
fid=fopen(fullfile(out,[name,'_metrics.json']),'w');fprintf(fid,'%s',jsonencode(r.metrics));fclose(fid);
u=r.velocity{1};v=r.velocity{2};speed=zeros(size(u));
for d=1:numel(r.velocity), speed=speed+r.velocity{d}.^2; end
speed=sqrt(speed);
if numel(r.edges)==3
    k=floor(size(speed,3)/2)+1;u=u(:,:,k);v=v(:,:,k);speed=speed(:,:,k);
end
[X,Y]=ndgrid(r.centers{1},r.centers{2});
f=figure('Visible','off');surf(X,Y,zeros(size(X)),speed,'EdgeColor','none');view(2);hold on;
quiver(X,Y,u,v,'w');axis equal tight;colorbar;xlabel('x (m)');ylabel('y (m)');title(strrep(name,'_',' '));
print(f,fullfile(out,[name,'.png']),'-dpng','-r150');close(f);
if numel(r.edges)==3
    f=figure('Visible','off');hold on;ex=r.edges{1};ey=r.edges{2};ez=r.edges{3};
    for z=[ez(1),ez(end)]
        for x=ex, plot3([x,x],[ey(1),ey(end)],[z,z],'b-');end
        for y=ey, plot3([ex(1),ex(end)],[y,y],[z,z],'b-');end
    end
    for y=[ey(1),ey(end)]
        for x=ex, plot3([x,x],[y,y],[ez(1),ez(end)],'b-');end
        for z=ez, plot3([ex(1),ex(end)],[y,y],[z,z],'b-');end
    end
    view(3);axis equal;xlabel('x');ylabel('y');zlabel('z');title('3D tensor mesh');
    print(f,fullfile(out,[name,'_mesh.png']),'-dpng','-r150');close(f);
end
end

function e=refine_grid(r,add)
dim=numel(r.edges);activity=zeros(size(r.pressure));
for comp=1:dim
    for d=1:dim
        g=derivative(r.velocity{comp},r.centers{d},d,dim);activity=activity+g.^2;
    end
end
e=cell(1,dim);
for d=1:dim
    weighted=activity.*r.volumes;
    for k=dim:-1:1, if k~=d, weighted=sum(weighted,k); end, end
    old=r.edges{d};score=weighted(:)'.*diff(old).^2;
    [~,order]=sort(score,'descend');chosen=order(1:min(add,numel(order)));
    e{d}=sort([old,(old(chosen)+old(chosen+1))/2]);
end
end

function g=derivative(u,x,d,dim)
% Three-point derivative on nonuniform coordinates, including end points.
order=[d,setdiff(1:dim,d,'stable')];v=permute(u,order);sz=size(v);a=reshape(v,numel(x),[]);b=zeros(size(a));
for i=1:numel(x)
    if i==1, ids=1:3;elseif i==numel(x),ids=numel(x)-2:numel(x);else,ids=i-1:i+1;end
    for j=1:3
        other=setdiff(1:3,j); w=(2*x(i)-x(ids(other(1)))-x(ids(other(2))))/((x(ids(j))-x(ids(other(1))))*(x(ids(j))-x(ids(other(2)))));
        b(i,:)=b(i,:)+w*a(ids(j),:);
    end
end
g=ipermute(reshape(b,sz),order);
end

function err=sample_error(r,reference)
dim=numel(r.edges);axes=cell(1,dim);query=axes;
for d=1:dim
    e=r.edges{d};axes{d}=linspace(e(1)+.15*(e(end)-e(1)),e(1)+.85*(e(end)-e(1)),9);
end
[query{:}]=ndgrid(axes{:});numerator=0;denominator=0;
for d=1:dim
    a=interpn(r.centers{:},r.velocity{d},query{:},'linear');
    b=interpn(reference.centers{:},reference.velocity{d},query{:},'linear');
    numerator=numerator+sum((a(:)-b(:)).^2);denominator=denominator+sum(b(:).^2);
end
err=sqrt(numerator/denominator);
end
