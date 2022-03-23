function dd = integrate_with_XSPDE(X0, F, G, opts)
% [F,G] = make_update_function_VO2_tdeptemp(opts);
% use the xSPDE solver 
% these are internal parameters for xSPDE
clear r
r.dimension= 2;
r.points = [ opts.npoints(1), 4*prod(opts.npoints(2:end)) ];
r.ranges = [ opts.ranges(1) prod(opts.ranges(2:end)) ];
r.noises  = 1;
r.fields  = 1;
r.seed    = randi(2^16);
r.da      = @(X, w, r) F(X, r.t)+G(X,r.t).*w;
r.initial = @(rv, w) X0;

% run `xsim` from xSDPE and extract the data from the structure `d`
[~,~,d] = xsim(r);

d = squeeze(d{1}{1});
dd = reshape(d(:,:,1), [opts.npoints 4]);

end

