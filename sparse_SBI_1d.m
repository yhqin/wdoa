function sp  = sparse_SBI_1d(R, n, design, wavelength, grid_size, lambda, varargin)
%Syntax:
%   sp = sparse_SBI_1d(R, n, design, wavelength, grid_size, lambda, ...);
%Inputs:
%   R - Sample covariance matrix.
%   design - Array design. Can also be a function handle that generates
%            a steering matrix. This function must take two arguments,
%            wavelength and the doa vector.
%   n - Number of sources.
%   wavelength - Wavelength. If design is set to a function handle, this
%                parameter must be set to [].
%   grid_size - Number of grid points used.
%   lambda - 
%Output:
%   sp - Spectrum structure with the following fields:
%           x - An 1 x grid_size vector.
%           y - An 1 x grid_size vector. Calling `plot(x, y)` will plot the
%               spectrum.
%           x_est - An 1 x n vector storing the estimated DOAs. May not
%                   fall on the grid if 'RefineEstimates' is set to true.
%           x_unit - The same as the unit specified by 'Unit'.


unit = 'radian';
is_noise_var_known = false;
noise_var = -1;
use_constrained_formulation = false;
upper_bound_l2_norm = false;
verbose = false;
for ii = 1:2:nargin-6
    option_name = varargin{ii};
    option_value = varargin{ii+1};
end

% discretize and create the corresponding steering matrix
resolution = 1
[doa_grid, doa_grid_display, ~] = default_doa_grid(grid_size, unit, resolution);
if ishandle(design)
    A = design(wavelength, doa_grid);
else
    A = steering_matrix(design, wavelength, doa_grid);
end
% preparing the objective function
[m, ~] = size(A);
if is_noise_var_known && ~upper_bound_l2_norm
    % remove noise variance from the main diagonal if known
    R = R - noise_var*eye(m);
end
r = R(:);
Phi = khatri_rao(conj(A), A);

 [xsize, ysize] = size(Phi);
 for ii = 1:xsize
     for jj = 1:ysize
        B(ii,jj) = -1i * pi * (ii-(design.element_count+1)/2) * sin(doa_grid(jj)/180*pi) * Phi(ii,jj);
     end
 end
 
  Phi_full = [real(Phi); imag(Phi)];
r_full = [real(r); imag(r)];
B_full = [real(B); imag(B)];
[dim_s, dim_x] = size(Phi_full);

 params.Y = r_full;
 params.A = Phi_full;
 params.B = B_full;
 params.resolution = resolution/180*pi;
 params.rho = 1e-2;
 params.alpha = mean(abs(Phi_full'*r_full), 2);
 params.beta = zeros(dim_x,1);
 params.K = n;
 params.maxiter = 2000;
 params.tolerance = 1e-3;
 params.sigma2 = mean(var(r_full))/100;


 el_time = cputime;
 res = SBI(params);  %SOMP_AltDesc(Phi_full,B_full,r_full,n,'SOMP-TLS',maxiters);
 el_time = cputime - el_time;
 x = res.mu;
 

if isempty(x)
    x = zeros(dim_x, 1);
end

% find n DOAs
[x_est_intl, ~, resolved] = find_doa_from_spectrum_1d(doa_grid_display, x, n);
% return
sp = struct();
sp.x = doa_grid_display;
sp.x_est = x_est_intl;
sp.x_unit = unit;
sp.y = x';
sp.resolved = resolved;
sp.discrete = true;
end

