function nllh = garchNllh(params, y)
%

sigma0 = 1;
k = params(1);
garch = params(2);
arch = params(3);

retrieveSigmas = zeros(numel(y), 1);
retrieveSigmas(1) = sigma0;

for ii=2:numel(y)
    retrieveSigmas(ii) = sqrt(k + garch*retrieveSigmas(ii-1).^2 +...
        arch*y(ii-1).^2);
end


nllh = sum(0.5*log(retrieveSigmas.^2*2*pi) + ...
    0.5*(y.^2./retrieveSigmas.^2));

end