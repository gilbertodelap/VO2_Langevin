function [F,G] = make_update_function_VO2_tdeptemp(opts)
% 
a = opts.a;
r = opts.r;
sig = opts.sig; 
gam = opts.gam;   % damping 
Laplacian = opts.Laplacian;


Ndims = length(opts.npoints)-1;
for ii=1:Ndims
    x_ = 0:(opts.ranges(ii+1)/opts.npoints(ii+1)):opts.ranges(ii+1); x_=x_(1:end-1);
    xs{ii} = x_(:);
end

function du = drift_function(u, t)
    u = reshape(u, [], 4);
    q1 = u(:,1);
    p1 = u(:,2);
    q2 = u(:,3);
    p2 = u(:,4);

    tmp1 = Laplacian(q1);
    tmp2 = Laplacian(q2);

    du(:,1) = a*p1;
    du(:,3) = a*p2;
    du(:,2) = ( -r(xs,t) - (q1.^2 ) ).*q1 + 1*tmp1 - gam*a*p1;
    du(:,4) = ( -r(xs,t) - (q2.^2 ) ).*q2 + 1*tmp2 - gam*a*p2;
    du = du(:)';
end

function B = noise_function(u, t)
    B = zeros(size(u,2)/4, 4);
    B(:,2) = sig(xs,t);
    B(:,4) = sig(xs,t);
    B = B(:)';
end

F = @drift_function;
G = @noise_function;

end 
