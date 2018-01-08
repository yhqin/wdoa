function sp = sparse_SBL_1d(R, n, design, wavelength, grid_size, lambda, varargin)
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
verbose = false;
for ii = 1:2:nargin-6
    option_name = varargin{ii};
    option_value = varargin{ii+1};
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

r = R(:);
Phi = khatri_rao(conj(A), A);

Phi_full = [real(Phi); imag(Phi)];
r_full = [real(r); imag(r)];
[dim_s, dim_x] = size(Phi_full);

    options = SBLSet();
    options.Nsource = n
    options.covergence.error = 10^(-4)
    tic;
    
    [x, report] = SBL_v3p12( Phi_full, r_full, options );
    toc;

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

