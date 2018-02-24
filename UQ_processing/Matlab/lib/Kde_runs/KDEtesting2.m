
function KDEtesting2

% Use Metropoli Hasting to generate samples

N = 1000; % number of movement

% target distribution
mu  = [2,2];
sig = [2,0.5
    0.5,2];
p = @(x) 1/(2*pi)/sqrt(det(sig))*exp(-.5*(x-mu)/(sig)*(x-mu)');

% proposal distribution
% it is normal with
% covariance is b^2*eye(2) where b is
b = 1.2;

% pre allocate memory
x_out = zeros(N,2);
rej_flag = zeros(N,1);
x_rej = zeros(size(x_out));

x_out(1,:) = randn(1,2)*b + mu; % start point
for i = 2:N
    y = randn(1,2)*b + x_out(i-1,:); % sample from q ~ N(x(n-1),b^2)
    alpha = min(p(y)/p(x_out(i-1,:)),1);
    u = rand; % sample u from U[0,1)
    if(u<=alpha) % accept
        x_out(i,:) = y;
        rej_flag(i) = 1;
    else % reject
        x_out(i,:) = x_out(i-1,:);
        x_rej(i,:) = y;
    end
end

figure(fid)
plot(x_out(:,1),x_out(:,2),'o')
fid=fid+1;

p = kde( x_out', [0.3;0.3]);

figure(fid)
plot(p);
fid = fid+1;

f=figure(fid)
pcolor(hist(p));
shading flat
fid = fid+1;
saveas(f,'tsts.pdf')

figure(fid)
mesh(hist(p));

b=marginal( p, 1)

