function sp = sparse_SOMP_LS_1d(R, n, design, wavelength, grid_size, lambda, varargin)
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
%           resolved - True if the number of peaks in the spectrum is
%                      greater or equal to the number of sources.
%           discrete - Constant value true.

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
if upper_bound_l2_norm 
    warning('Specified noise variance will be ignore when the problem formulation is set to ''ConstrainedL2''.');
end
% discretize and create the corresponding steering matrix
[doa_grid, doa_grid_display, ~] = default_doa_grid(grid_size, unit, 1);
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

issvd = 1;
 if issvd
      [s,v,d] = svd(r,'econ');
       r = r * d    % d(:, 1:design.element_count);
 end
  % Perform SOMP-LS Alternate Descent
  %-------------------------------------------------------------%
 maxiters = 2*n;
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
 el_time = cputime;
 [D1,x] = SOMP_AltDesc(Phi_full,B_full,r_full,n,'SOMP-LS',maxiters);
 el_time = cputime - el_time;
            
 [~,supp1] = sort(sum(abs(x).^2,2),'descend');supp1 = supp1(1:n)';
 D1 = D1(supp1,supp1);

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

