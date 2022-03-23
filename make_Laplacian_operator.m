function Laplacian = make_Laplacian_operator(opts)

clengths2 = opts.cohlengths2;

Ndims = length(opts.npoints)-1;

% the first coordinate (x) is OBC, the rest at Periodic BC
ks{1} = (0:(opts.npoints(2)-1))*pi/opts.ranges(2);
for ii=2:Ndims
    ks{ii} = (-opts.npoints(ii+1)/2:(opts.npoints(ii+1)/2-1))*2*pi/opts.ranges(ii+1);
    ks{ii} = fftshift(ks{ii});
%     N = opts.npoints(ii+1);
%     L = opts.ranges(ii+1);
%     ks{11} = (-N/2:(N/2-1))*2*pi/L;
%     ks{ii} = fftshift(ks{ii});
end

if Ndims == 1
    epsilon2 = clengths2(1)*ks{1}.^2;
elseif Ndims == 2
    [kx,ky] = ndgrid(ks{1}, ks{2});
    epsilon2 = (clengths2(1)*kx.^2 + clengths2(2)*ky.^2);
elseif Ndims == 3
    [kx,ky,kz] = ndgrid(ks{1}, ks{2}, ks{3});
    epsilon2 = clengths2(1)*kx.^2 + clengths2(2)*ky.^2 + clengths2(3)*kz.^2;
%     epsilon2 = kx.^2 + ky.^2 + kz.^2;
end
epsilon2 = -epsilon2(:);

function uk = bc(u)
    u = reshape(u, [opts.npoints(2:end), 1]);
    uk = dct(u,[], 1);
    for d = 2:Ndims
        uk = fft(uk,[],d);
    end
    uk = uk(:);
end

function u = ibc(uk)
    uk = reshape(uk, [opts.npoints(2:end), 1]);
    u = idct(uk, [], 1);
    for d = 2:Ndims
        u = ifft(u,[],d);
    end
    u = u(:);
end

function u1 = Laplacian_operator(u)
    u1 = ibc(epsilon2.*bc(u));
end

Laplacian = @Laplacian_operator;

end
