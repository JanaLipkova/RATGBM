function M1 = resize3d(M,new_size,interp)
%resize 3D images

% check for interpolation mode input
if nargin == 2
    interp_mode = '*linear';
elseif nargin ~= 3
    error('incorrect number of inputs');
else
    interp_mode = interp;
end



% extract the sizes
[Nx_input, Ny_input, Nz_input] = size(M);
Nx_output=new_size(1);
Ny_output=new_size(2);
Nz_output=new_size(3);
       

% update command line status
disp(['  input grid size: ' num2str(Nx_input) ' by ' num2str(Ny_input) ' by ' num2str(Nz_input) ' elements']);
disp(['  output grid size: ' num2str(Nx_output) ' by ' num2str(Ny_output) ' by ' num2str(Nz_output) ' elements']); 

% create normalised plaid grids of current discretisation
[x_mat, y_mat, z_mat] = ndgrid((0:Nx_input-1)/(Nx_input-1), (0:Ny_input-1)/(Ny_input-1), (0:Nz_input-1)/(Nz_input-1));       

% create plaid grids of desired discretisation
[x_mat_interp, y_mat_interp, z_mat_interp] = ndgrid((0:Nx_output-1)/(Nx_output-1), (0:Ny_output-1)/(Ny_output-1), (0:Nz_output-1)/(Nz_output-1));

% compute interpolation; for a matrix indexed as [M, N, P], the
% axis variables must be given in the order N, M, P
M1 = interp3(y_mat, x_mat, z_mat, M, y_mat_interp, x_mat_interp, z_mat_interp, interp_mode);        

