function run_case()
% Periodic conservative advection-diffusion; same case as Python.
base=fileparts(fileparts(mfilename('fullpath')));
c=jsondecode(fileread(fullfile(base,'cases','transport.json')));
out=fullfile(base,'results','matlab'); if ~exist(out,'dir'), mkdir(out); end
schemes={'upwind','lax_wendroff'}; rows=struct([]); q=0;
f=figure('Visible','off'); hold on;
for s=1:2
    for ni=1:numel(c.cell_counts)
        n=c.cell_counts(ni); dx=c.length/n; x=((0:n-1)'+.5)*dx;
        assert(n>=8 && c.diffusivity>=0 && c.final_time>0 && c.length>0);
        rate=abs(c.speed)/dx+2*c.diffusivity/dx^2;
        steps=max(1,ceil(c.final_time*rate/.8)); dt=c.final_time/steps;
        C=c.speed*dt/dx; D=c.diffusivity*dt/dx^2;
        assert(C^2+2*D<=1+1e-12);
        u=1+.5*sin(2*pi*x/c.length); mass=dx*sum(u); timer=tic;
        for j=1:steps
            right=circshift(u,-1);
            if s==1
                if c.speed>=0, upstream=u; else, upstream=right; end
                flux=c.speed*upstream-c.diffusivity*(right-u)/dx;
            else
                flux=c.speed*.5*(u+right)-.5*c.speed*C*(right-u)-c.diffusivity*(right-u)/dx;
            end
            u=u-dt/dx*(flux-circshift(flux,1));
        end
        elapsed=toc(timer);
        exact=1+.5*exp(-c.diffusivity*(2*pi/c.length)^2*c.final_time)*sin(2*pi*(x-c.speed*c.final_time)/c.length);
        q=q+1; rows(q)=struct('cells',n,'scheme',schemes{s},'dt',dt,'steps',steps,'courant',C,'diffusion_number',D,...
            'l2_error',sqrt(dx*sum((u-exact).^2)),'mass_drift',dx*sum(u)-mass,'integration_seconds',elapsed);
        writetable(table(x,u,exact,'VariableNames',{'x','numerical','analytical'}),fullfile(out,sprintf('%s_%d.csv',schemes{s},n)));
    end
    plot(x,u,'DisplayName',schemes{s});
end
plot(x,exact,'k--','DisplayName','Analytical');legend('show','Interpreter','none');
xlabel('x (m)');ylabel('Scalar');title('Periodic advection-diffusion');
print(f,fullfile(out,'solution.png'),'-dpng','-r150');close(f);
writetable(struct2table(rows),fullfile(out,'convergence.csv'));save(fullfile(out,'run.mat'),'c','rows');disp(struct2table(rows));
end
