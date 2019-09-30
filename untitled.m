i = 0;
% for j = 70:150
Currents1 = z.grid{1,3}(84,:);
Currents2 = z.grid{1,3}(85,:);
kk = 0.85;
Currents3 = kk*Currents1 + (1-kk)*Currents2;
plist.Currents = Currents3;
[psi,Dpsi,x_point,m_axis,ax_val,bdry_val,id_bnd,id_ax,Meshes,varargout] = ...
    solve_FreeBoundaryAdaptive(plist,Mesh,AdaptDepth,Output_type);
% [psi,Dpsi,x_point,m_axis,ax_val,bdry_val,id_bnd,id_ax] = ...
%     solve_FreeBoundary(plist,Mesh);
% if varargout ~=0
%     i=i+1;
%     marker(i) = j;
% end
% end

psi_idbnd = psi(id_bnd);
psi = project2grid(psi,Meshes,CommonMesh);
output = [psi;x_point';m_axis';psi_idbnd];

contour_plot(psi,CommonMesh,1,psi_idbnd,'red','--',0.5);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Currents CommonMesh
plist = MyReadYaml('~/ITER_IPfree.yml');
Currents = plist.Currents; % Store current values
mat_file = '../Meshdata/ITER.mesh.mat'; %[plist.mesh_file, '.mat'];
Mesh = load(mat_file);
CommonMesh = RefineInterior(Mesh,2); % Common Mesh to project the adaptive solution
QoI_l = length(CommonMesh.Coordinates)+5;
noise_level = 0.01;  % Noise level
ColDepth = 2;  % Maximum number of levels to compute collocation nodes
AdaptDepth = 1; 
Output_type = 'LastLevOut';

%%%---------------------------
figure
T = importdata('fail_count.txt');
index = find(T~=0) - 25;
% BadCurrents = zeros(length(index), 12);
for i = 1:length(index)
%     BadCurrents(i,:) = z.grid{3}(index(i),:);
    
%     plist.Currents = BadCurrents(i,:);
%     [psi,Dpsi,x_point,m_axis,ax_val,bdry_val,id_bnd,id_ax,Meshes,varargout] = ...
%     solve_FreeBoundaryAdaptive(plist,Mesh,AdaptDepth,Output_type);
% psi_idbnd = psi(id_bnd);
% psi = project2grid(psi,Meshes,CommonMesh);
psi = z.fvals{3}(index(i),1:QoI_l - 5);
psi_idbnd = z.fvals{3}(index(i),end);
            contour_plot(psi',CommonMesh,1,psi_idbnd,'red','--',0.5);

end
