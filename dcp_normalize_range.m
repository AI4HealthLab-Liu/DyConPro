function normdat=dcp_normalize_range(X,rlow,rhigh)

% Normalize range to interval [0 1], or to interval [rlow rhigh]


minx=min(X);
maxx=max(X);
rang=maxx-minx;

normdat=(X-minx)./rang; % norm data to interval [0 1]; there will be a 0 and a 1 in each vector

% norm data to range interval [rlow rhigh]; there will be a rlow and rhigh value in each vector
if ~isempty(rlow) && ~isempty(rhigh)
    rang2=rhigh-rlow;
    normdat=(normdat.*rang2)+rlow;
end



end